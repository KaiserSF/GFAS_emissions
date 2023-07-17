import json
import gribapi as ga
import pandas as pd

# dominant land cover class according to map values -> [0:7]
LCC = ('SA', 'SAOS', 'AG', 'AGOS', 'TF', 'PEAT', 'EF', 'EFOS')

def read_ef_cvs(fname, write_json=False):
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
    spec_report_to_grib = {
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
        'Formaldehyde'  : 'ch2oh'    ,
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
    #for k in spec_report_to_grib.keys():
    #    spec_report_to_grib[k] = spec_report_to_grib[k] + 'fire'
    df.insert(len(df.columns), 'shortName', [spec_report_to_grib[x] for x in df['Species']])

    print('building dictionary all emission factors...')
    ef = {}
    for v in ('v1p2','a19'):
        ef[v] = {}
        for s in spec_report_to_grib.values():
            ef[v][s] = {}
            for l in lc_a19_to_v1p2.values():
                #print(v,l,s)
                #print (df.loc[(df['lc']==l) & (df['shortName']==s)])
                ef[v][s][l] = df.loc[(df['lc']==l) & (df['shortName']==s)].iloc[0][v]

        print(f'adding land cover types without distince emission factors, i.e. SAOS, AGOS got {v}...')
        for s in spec_report_to_grib.values():
            ef[v][s]['SAOS'] = ef[v][s]['SA']
            ef[v][s]['AGOS'] = ef[v][s]['AG']

    if write_json:
        json_file = 'emission_factors.json'
        print(f'writing emission factors to {json_file} ...')
        with open(json_file, 'w', encoding='utf-8') as fp:
            json.dump(ef, fp, ensure_ascii=False, indent=2)
        
    return ef
    
def GFAS_emissions(dm_fname, emi_fname, dlc_fname, ef_fname):
    
    # read field of combustion rate, resp. dry matter burnt
    with open(dm_fname, 'r') as fp:
        dm = ga.grib_new_from_file(fp)
    assert(ga.grib_get(dm, 'shortName') == 'crfire')
    print(ga.grib_get(dm, 'paramId'), ga.grib_get(dm, 'shortName'), ga.grib_get(dm, 'name'))

    # read land cover map
    with open(dlc_fname, 'r') as fp:
        dlc = ga.grib_new_from_file(fp)
    assert(ga.grib_get(dlc, 'paramId') == 94)

    # read species emission factors
    ef = read_ef_cvs(ef_fname)

    # empty output file
    with open(emi_fname, 'w') as fp:
        pass

    # loop over species
    for s in ['nox','co','nh3']:
        
        # calculate species emission fluxes
        emi = ga.grib_clone(dm)
        ga.grib_set(emi, 'shortName', s+'fire')
        print(ga.grib_get(emi, 'paramId'), ga.grib_get(emi, 'shortName'), ga.grib_get(emi, 'name'))

        # write species emission fluxes
        with open(emi_fname, 'ab') as fp:
            ga.grib_write(emi, fp)
    
if __name__ == '__main__':
    GFAS_emissions('fields.grb', 'emi.grb', 'dat/dlc.grb', 'dat/Table2_GFAS_vs_A19_EF_summary_longformat.csv')
