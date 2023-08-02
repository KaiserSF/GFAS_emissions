import json
import gribapi as ga
import pandas as pd
import numpy as np

# dominant land cover class according to map -> map index to [1:8]!
LCC = ('SA', 'SAOS', 'AG', 'AGOS', 'TF', 'PEAT', 'EF', 'EFOS')
SPEC_REPORT_TO_GRIB = {
        'Acetaldehyde'  : 'c2h4o'    ,
        'Acetone'       : 'c3h6o'    ,
        'BC_or_EC'      : 'bc'       ,
        'Benzene'       : 'c6h6'     ,
        'Butanes'       : 'c4h10'    ,
        'Butenes'       : 'c4h8'     ,
        'C'             : 'c'        ,
        'C2H4'          : 'c2h4'     ,
        'C2H6'          : 'c2h6'     ,
        'C3H6'          : 'c3h6'     ,
        'C3H8'          : 'c3h8'     ,
        'CH4'           : 'ch4'      ,
        'CO'            : 'co'       ,
        'CO2'           : 'co2'      ,
        'DMS'           : 'c2h6s'    ,
        'Ethanol'       : 'c2h5oh'   ,
        'Formaldehyde'  : 'ch2o'    ,
        'H2'            : 'h2'       ,
        'Heptanes'      : 'c7h16'    ,
        'Hexanes'       : 'c6h14'    ,
        'Hexene'        : 'c6h12'    ,
        'Higher_Alkanes': 'hialkanes',
        'Higher_Alkenes': 'hialkenes',
        'Isoprene'      : 'c5h8'     ,
        'Methanol'      : 'ch3oh'    ,
        'N2O'           : 'n2o'      ,
        'NH3'           : 'nh3'      ,
        'NMHC_sum'      : 'nmhc'     ,
        'NOx_as_NO'     : 'nox'      ,
        'OC'            : 'oc'       ,
        'Octenes'       : 'c8h16'    ,
        'PM2.5'         : 'pm2p5'    ,
        'Pentane'       : 'c5h12'    ,
        'Pentenes'      : 'c5h10'    ,
        'SO2'           : 'so2'      ,
        'TC'            : 'tc'       ,
        'TPM'           : 'tpm'      ,
        'Terpenes'      : 'terpenes' ,
        'Toluene'       : 'c7h8'     ,
        'Toluene_lump'  : 'toluene'  ,
        'Xylenes'       : 'c8h10'
        }

def read_ef_cvs(fname, write_json=True):
    df = pd.read_csv(fname)
    print(f'read {fname} with {df.columns.values}')
    right_cols = (df.columns.values==['Species','A19_average','gfas_v1p2','perc_change','type'])
    assert(right_cols.all())

    df = df.rename(columns={'A19_average':'a19','gfas_v1p2':'v1p2'})
    print(f'new column names: {df.columns.values}')

    print('mapping land cover classes of Andreae 2019 to GFASv1.2: temperate forest is EF, boreal forest is EFOS!')
    lc_a19_to_v1p2 = {
        'peat'                  : 'PEAT',
        'tropical_forest'       : 'TF',
        'temperate_forest'      : 'EF',
        'savannah_and_grassland': 'SA',
        'boreal_forest'         : 'EFOS',
        'agricultural_residue'  : 'AG' }
    df.insert(len(df.columns), 'lc', [lc_a19_to_v1p2[x] for x in df['type']])

    print('mapping species names to GRIB short names...')
    df.insert(len(df.columns), 'shortName', [SPEC_REPORT_TO_GRIB[x] for x in df['Species']])

    print('building dictionary all emission factors...')
    ef = {}
    for v in ('v1p2','a19'):
        ef[v] = {}
        for s in SPEC_REPORT_TO_GRIB.values():
            ef[v][s] = {}
            for l in lc_a19_to_v1p2.values():
                #print(v,l,s)
                #print (df.loc[(df['lc']==l) & (df['shortName']==s)])
                ef[v][s][l] = df.loc[(df['lc']==l) & (df['shortName']==s)].iloc[0][v] / 1000

        print(f'adding land cover types without distince emission factors, i.e. SAOS, AGOS got {v}...')
        for s in SPEC_REPORT_TO_GRIB.values():
            ef[v][s]['SAOS'] = ef[v][s]['SA']
            ef[v][s]['AGOS'] = ef[v][s]['AG']

    if write_json:
        json_file = 'emission_factors.json'
        print(f'writing emission factors to {json_file} ...')
        with open(json_file, 'w', encoding='utf-8') as fp:
            json.dump(ef, fp, ensure_ascii=False, indent=2)
        
    return ef
    
def GFAS_emissions(dm_fname, dlc_fname, ef_fname, action='emission'):
    '''
    calculate various species emissions from dry matter burnt field (GRIB shortName="crfire")
    based on emission factors from Andreae & Merlet 2001 and Andreae 2019
    and/or calculate ratios of the emission factors (A&M2019 / A2001)
    
    Arguments:
        dm_fname: file name of grib input field of dry matter burnt
        dlc_fname: file name of GFAS land cover field
        ef_name: file name of emission factor table by CAMS_44 (de Jong et al.)
        action: "emission" or "ratio" or combination therefore to determine calculation
    Output:
        emissions and/or ratios for each land cover type and chemical species in various files
    Return:
        None
    '''
    # read land cover map
    with open(dlc_fname, 'r') as fp:
        dlcmsg = ga.grib_new_from_file(fp)
    assert(ga.grib_get(dlcmsg, 'paramId') == 94)
    dlc = ga.grib_get_values(dlcmsg)

    # read species emission factors
    ef = read_ef_cvs(ef_fname)

    if 'ratio' in action:
        # RATIOS OF EFS
        
        # loop over land cover types
        for il, l in enumerate(LCC):
        
            # empty output file
            out_fname = f'Andreae2019_over_GFASv1p2_{l}.grb'
            with open(out_fname, 'w') as fp:
                pass
        
            # loop over species
            for s in ['nox','co','nh3']:
            
                # calculate ratios of emission factors
                emi = ga.grib_clone(dlcmsg)
                ga.grib_set(emi, 'shortName', s+'fire')
                #print(f'processing {ga.grib_get(emi,"paramId")} - {ga.grib_get(emi,"shortName")} in {l}...')
                val = ga.grib_get_values(emi)
                ratio = ef['a19'][s][l] / ef['v1p2'][s][l]
                print(f'ratio of {s} in {l}: {ratio}')
                val = np.where(np.abs(il+1-dlc)<.1, ratio, 1)
                ga.grib_set_values(emi, val)
                
                # write species emission fluxes
                with open(out_fname, 'ab') as fp:
                    ga.grib_write(emi, fp)

    if 'emission' in action:
        # EMISSIONS
        
        # read field of combustion rate, resp. dry matter burnt
        with open(dm_fname, 'r') as fp:
            dmmsg = ga.grib_new_from_file(fp)
        assert(ga.grib_get(dmmsg, 'shortName') == 'crfire')
        print(ga.grib_get(dmmsg, 'paramId'), ga.grib_get(dmmsg, 'shortName'), ga.grib_get(dmmsg, 'name'))
        dm = ga.grib_get_values(dmmsg)
        
        # loop over emission factor set versions
        for v in ('a19','v1p2'):
        
            # empty output file
            out_fname = f'emissions_GFAS_{v}.grb'
            with open(out_fname, 'w') as fp:
                pass
        
            # loop over species
            for s in SPEC_REPORT_TO_GRIB.values():
                print(f'processing {v},{s}')
                
                # calculate ratios of emission factors
                emi = ga.grib_clone(dlcmsg)
                ga.grib_set(emi, 'shortName', s+'fire')
                val = ga.grib_get_values(emi)
                #print(val)
                val[:] = 0
                #print(val)
                  
                # loop over land cover types
                for il, l in enumerate(LCC):
                    val += np.where(np.abs(il+1-dlc)<.1, dm*ef[v][s][l], 0)
   
                ga.grib_set_values(emi, val)
               
                # write species emission fluxes
                with open(out_fname, 'ab') as fp:
                    ga.grib_write(emi, fp)
                    
    return None

if __name__ == '__main__':
    # The first argument is the input field of dry matter burnt, resp. crfire.
    GFAS_emissions('ADSfields.grib', 'dat/dlc.grb', 'dat/Table2_GFAS_vs_A19_EF_summary_longformat.csv', 'ratio+emissions')
