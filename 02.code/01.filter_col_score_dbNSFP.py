# 2025.01.19
# by goslak
# 01.filter_col_score_dbNSFP.py

# canonical : most conserved, most highly expressed, the longest coding sequence

# import library ##############################################################################
import pandas as pd

# read_file ##############################################################################
def read_file(_para):
    print("===[ Def: read_file ]===")
    _in = para['inFold'] + '/' + para['inFile']
    _tmp = pd.read_csv(_in, sep="\t", index_col = False, low_memory=False)
    print("!!! Read :", _in)
    print("\tN =", _tmp.shape)
    print(_tmp.head(1))
    return(_tmp)

_prefix = 'clinvar2021'
_prefix = 'clinvar2017_2020.47394'

para = {'inFold': '/home/dna/99.github/compare_pathogenic_tool/00.data',
        'inFile': f'{_prefix}.dbNSFP.out'}
df = read_file(para)

para['inFile'] = f'{_prefix}.dataset'
dataset = read_file(para)
dataset.rename(columns={'pos': 'pos(1-based)', 'alRef': 'ref', 'alAlt': 'alt'}, inplace=True)
print(dataset.head(1))

# filter_column ##############################################################################
def filter_column(_df, _para):
    print("===[ Def: filter_column ]===")
    
    _tmp = _df.iloc[:, _para['keep_col']].copy()
    _tmp['flag'] = 0
    for _v, _cols in _para['check_v'].items():
        for _col in _cols:
            _count = 0
            if _v == '.':
                _tmp.loc[_tmp[_col] == _v, 'flag'] = 1
                _count = (_tmp[_col] == _v).sum()
            elif _v == 'YES':
                _count = (~_tmp[_col].str.contains(_v, na=False)).sum()
                _tmp.loc[~_tmp[_col].str.contains(_v, na=False), 'flag'] = 1
            else:
                pass
            print(f'\t{_col}\t{_v}\t{_count}')
    print(f"\tDrop(flag=1)\t{(_tmp['flag'] == 1).sum()}")
    
    _tmp = _tmp[_tmp['flag'] != 1].reset_index(drop=True)
    _tmp = _tmp.drop('flag', axis=1)
    
    _cols = [_col for _col in _tmp.columns 
                if _col.endswith(('rankscore','group','protID','AAchange','Top5features','model','AAE','LOO','0_pred')) or 
                _col.startswith(('Aloft','GM12878','H1','HUVEC','integrated')) or
                #'phred' in _col or 'LOO' in _col or 'Class25' in _col or '-PC' in _col]
                'phred' in _col or 'Class25' in _col]
    _tmp = _tmp.drop(_cols, axis=1)
    print("\tN =", _tmp.shape)
    print(_tmp.head(1))
    return(_tmp)

para['keep_col'] = (list(range(0,6)) + list(range(11,16)) + list(range(27,31)) + 
                    list(range(460,465)) + [439,440,442] +
                    list(range(37,75)) + list(range(76,155)) + list(range(158,211)) + 
                    [212,224,226,228,230,232,234,284,361])
para['check_v'] = {'.':('aaref','aaalt'), 'YES':('VEP_canonical',)}
df_f = filter_column(df, para)

# save_file ##############################################################################
print("===[ save_file ]===")

para['outFile'] = f'{_prefix}.dbNSFP.out.filter'
_out = para['inFold'] + '/' + para['outFile']

print("!!! Save : ", _out)
print("\tN =", df_f.shape)
df_f.to_csv(_out, sep='\t', index=False,)

# check_canonical ##############################################################################
def check_canonical(_df, _col):
    print("===[ Def: check_canonical ]===")
    
    _df['yes'] = _df[_col].apply(lambda x: ','.join(str(i) for i, value in enumerate(x.split(';')) if value == 'YES'))
    print(_df['yes'].value_counts())
    return(_df)

df_f = check_canonical(df_f, 'VEP_canonical')

# eq_num ##############################################################################
def eq_num(_df, _para):
    print("===[ Def: eq_num ]===")
    _tmp = pd.DataFrame()
    _cols = _para['eqCols']

    for _idx, _col in enumerate(_cols):
        _tmp[_col] = _df[_col].apply(lambda x: len(x.split(';')))

        _tmp['eq'] = (_tmp[_cols[0]] == _tmp[_col])
        if _tmp['eq'].all():
            pass
        else:
            print(f"\tNo(eq) : {_cols[0]} & {_col}\t{_tmp[_tmp['eq']==False].shape}")
            print(_tmp[_tmp['eq']==False].head())
    
    _df['num'] = _tmp[_cols[0]]
    return(_df)

para['eqCols'] = (df_f.loc[:, 'aapos':'codonpos'].columns.tolist() + 
                df_f.loc[:, 'SIFT_score':'Polyphen2_HVAR_pred'].columns.tolist() +
                df_f.loc[:, 'MutationAssessor_score':'VEST4_score'].columns.tolist() + 
                ['MetaRNN_score','MetaRNN_pred'] + 
                df_f.loc[:, 'REVEL_score':'MPC_score'].columns.tolist() + 
                ['DEOGEN2_score','DEOGEN2_pred'] + 
                df_f.loc[:, 'LIST-S2_score':'PHACTboost_score'].columns.tolist())
#print(para['eqCols'])
df_f = eq_num(df_f, para)

# extract_score ##############################################################################
def extract_score(_df, _para):
    print("===[ Def: extract_score ]===")
    for _col in _para['eqCols']:
        _df[_col] = _df.apply(
            lambda row: row[_col].split(';')[int(row['yes'].split(',')[0])]
            if isinstance(row['yes'], str) else None,
            axis=1
        )
    return(_df)

df_s = extract_score(df_f, para)

# check_clinvar ##############################################################################
def check_clinvar(_df, _para):
    print("===[ Def: check_clinvar ]===")
    _clnsig = _para['clnsig']

    clnsig = {
        'Benign': 1,'Benign/Likely_benign': 2,'Likely_benign': 3,
        'Benign|association': 1,
        'Pathogenic': 4, 'Pathogenic/Likely_pathogenic': 5, 'Likely_pathogenic': 6, 
        'Uncertain_significance':-2,
        'Conflicting_classifications_of_pathogenicity':-1,
        #'drug_response': -1,
        #'Conflicting_classifications_of_pathogenicity|risk_factor': -1,
        #'not_provided': -2, '.': -2,
        #'Uncertain_significance|drug_response': -2,
        #'Likely_benign|association': 3,
        #'Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance': 5,
        #'Benign|drug_response': 1,
        #'Pathogenic/Likely_pathogenic/Likely_risk_allele': 5,
        #'Pathogenic/Likely_pathogenic|risk_factor': 5,
        #'Benign|other': 1,
        #'Pathogenic|risk_factor': 4,
        #'Likely_pathogenic/Likely_risk_allele': 6,
        #'Pathogenic|drug_response': 4,
        #'Likely_pathogenic|drug_response': 6,
        #'Benign/Likely_benign|other|risk_factor': 3,
        #'Benign/Likely_benign|association': 3,
        #'Benign/Likely_benign|other': 3,
        #'Pathogenic/Likely_risk_allele': 4
    }

    review = {
        'criteria_provided,_multiple_submitters,_no_conflicts': 'M',
        'criteria_provided,_conflicting_classifications': 'C',
        'criteria_provided,_single_submitter': 'S',
        'reviewed_by_expert_panel': 'E',
        'practice_guideline': 'G', #'.': 'U'
    }

    _df['clinvar_clnsig'] = _df['clinvar_clnsig'].map(clnsig).fillna(-3)
    _df['clinvar_review'] = _df['clinvar_review'].map(review).fillna('Unknown')
    
    _df[_clnsig] = 0
    _df.loc[(0 < _df['clinvar_clnsig']) & (_df['clinvar_clnsig'] <= 3), _clnsig] = 0
    _df.loc[(4 <= _df['clinvar_clnsig']) & (_df['clinvar_clnsig'] <= 6), _clnsig] = 1
    _df.loc[(0 > _df['clinvar_clnsig']), _clnsig] = -1

    print(_df['clinvar_clnsig'].value_counts().sort_index())
    print(_df['clinvar_review'].value_counts().sort_index())
    print(_df[_clnsig].value_counts().sort_index())
    return(_df)

para['clnsig'] = 'clnSig24'
df_s = check_clinvar(df_s, para)

# add varType ##############################################################################
def add_varType(_df, _para):
    print("===[ Def: add_varType ]===")
    _varType = _para['varType']
    _df['aapos'] = _df['aapos'].astype('int')
    _df[_varType]='m'
    _df.loc[_df['aapos'] == 1, _varType]='i'
    _df.loc[_df['aaref'] == 'X', _varType]='t'
    _df.loc[_df['aaalt'] == 'X', _varType]='n'
    
    print(_df[_varType].value_counts().sort_index())
    return(_df)

para['varType'] = 'varType24'
df_s = add_varType(df_s, para)

# check_dup ##############################################################################
def check_dup(_df, _para):
    print("===[ Def: check_dup ]===")
    print("### Data :",_para['data'])
    _dupIdx = ''
    for _k, _v in _para['dup'].items():
        print("!!! DUP : ",_v)
        _n = _df.duplicated(subset=_v).sum()
        print("\tN(dup) =", _n)
        if _n > 0:
            _tmp = _df[_df.duplicated(subset=_v, keep=False)]
            print("\tN(Row) =", _tmp.shape)
            if _k == 'Dup0' and _para['data'] == 'score':
                _dupIdx = _tmp[_tmp['CCDS_id']=='.'].index
                print("\tIndex(dup) =", _dupIdx)
            print(_tmp[_para['viewCols']].sort_values(by=['#chr','pos(1-based)']).head())
    
    return(_dupIdx)
  
para['dup'] = {'Dup0':['#chr', 'pos(1-based)','ref','alt'], 'Dup1':['#chr', 'pos(1-based)','ref']}
para['data'] = 'score'
para['viewCols'] = ['#chr','pos(1-based)','ref','alt', 'aaref','aaalt', 'aapos', 'genename', 'CCDS_id',
                    'clinvar_id','clinvar_clnsig',	'clinvar_review','varType24']
para['dupIdx'] = check_dup(df_s, para)
print("\tDrop(dup) =", len(para['dupIdx']))

if len(para['dupIdx']) > 0:
    df_s = df_s.drop(para['dupIdx']).reset_index(drop=True)
print("\tN =", df_s.shape)

para['data'] = 'dataset'
para['viewCols'] = ['#chr','pos(1-based)','strand','ref','alt','clnSig','clnStat','varType']
para['dupIdx'] = check_dup(dataset, para)         

# merge_file ##############################################################################
def merge_file(_df1, _df2, _para):
    print("===[ Def: merge_file ]===")
    _dupCols = _para['cols']
    _copyCols = _para['copyCols']
    _clnsig = _para['clnsig']
    _df1 = _df1.merge(_df2[_dupCols + _copyCols], on=_dupCols, how='left')

    _df1[_clnsig] = 0
    _df1.loc[_df1['2023_rs']=='B', _clnsig] = 0
    _df1.loc[_df1['2023_rs']=='P', _clnsig] = 1
    _df1.loc[~_df1['2023_rs'].isin(['B','P']), _clnsig] = -1

    print("\tN =",_df1.shape)
    print(_df1.head(1))
    print(pd.crosstab(_df1[_clnsig], _df1['2023_rs'], margins=True, margins_name="Total"))
    return(_df1)

para['cols'] = ['#chr', 'pos(1-based)', 'ref', 'alt']
para['copyCols'] = dataset.loc[:, 'sYear':'2023_rs'].columns.tolist() + ['varType']
para['clnsig'] = 'clnSig23'
df_s = merge_file(df_s, dataset, para)
#print(df_s.columns.tolist())

print("\tDiff(varType24 & varType) =",df_s[df_s['varType24']!=df_s['varType']].shape[0])
print(pd.crosstab(df_s['varType24'], df_s['varType'], margins=True, margins_name="Total"))

if df_s[df_s['varType24']!=df_s['varType']].shape[0] > 0:
    print(df_s[df_s['varType24']!=df_s['varType']][['#chr','pos(1-based)','varType24', 'varType']].head())
    df_s = df_s[df_s['varType24'] == df_s['varType']].reset_index(drop=True)
    print("\tN =", df_s.shape)

df_s = df_s.drop('varType24', axis=1)

# save_file ##############################################################################
print("===[ save_file ]===")

para['outFile'] = f'{_prefix}.dbNSFP.out.score'
_out = para['inFold'] + '/' + para['outFile']
print("!!! Save : ", _out)
print("\tN =", df_s.shape)

df_s = df_s.replace('.', '')
df_s.to_csv(_out, sep='\t', index=False,)

#print(df_s.columns.tolist())
print(pd.crosstab(df_s['clnSig24'], df_s['2023_rs'], margins=True, margins_name="Total"))
print(pd.crosstab(df_s['clnSig24'], df_s['clnSig23'], margins=True, margins_name="Total"))
