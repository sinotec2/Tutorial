from pandas import *
def unlist(InList):
    a=[]
    for i in InList:
        for j in i:
            a.append(j)
    return a
df=read_csv('ParamVals.csv')
df_ChNm=read_csv('ChemNames.csv')
sp_nam={'Antimony', 'Arsenic', 'Beryllium', 'Cadmium', 'Chromium(VI)', 'Cobalt',\
 'HeptaCDD, 1,2,3,4,6,7,8-', 'HexaCDD, 1,2,3,4,7,8-', 'HexaCDD, 1,2,3,6,7,8-',\
 'HexaCDD, 1,2,3,7,8,9 -', 'HexaCDD, 2,3,7,8-;', 'Lead', 'Manganese', 'Mercury','Methyl Mercury',\
 'Nickel', 'OctaCDD, 1,2,3,4,6,7,8,9-', 'PentaCDD, 1,2,3,7,8-', 'PentaCDD, 2,3,7,8-',\
 'Selenium', 'TetraCDD, 2,3,7,8-'}
sp_nam=list(sp_nam)
sp_cas=unlist([list(df_ChNm[df_ChNm['Name']==x]['CAS#']) for x in sp_nam])
#!grep Kow Definitions.csv
#69,"Kow (log)","unitless","Octanol-water partition coefficient","PC"
nam_K=['Kow','Kdsw','Kdbs','BSAFfish','BAF',"BCFfish",'Da',"Dw",'H','Kds']
i_K=[69,67,65,19,7,13,36,38,60,66]
d={}
df_final=DataFrame({'spnam':Series(sp_nam),'CAS':Series(sp_cas)})
for iK in i_K:
    namK=nam_K[i_K.index(iK)]
    df_K=df[df['Param']==iK]
    sp_boo=[x in set(df_K['CAS']) for x in sp_cas]
    df_spK=DataFrame({'spnam':Series(sp_nam),'CAS':Series(sp_cas),'withK':Series(sp_boo)})
    df_spK_Drop=df_spK[df_spK['withK']==True].reset_index()
    del df_spK_Drop['index']
    sp_K=unlist([list(df_K[df_K ['CAS']==x]['ParamVal']) for x in list(df_spK_Drop['CAS'])])
    df_spK_Drop[namK]=Series(sp_K)
    if namK=='Kow':
#http://publications.gc.ca/Collection/En1-34-2-2002E.pdf
#log(Kow) for Methy Mercury is between 1.72.5, take 2.1
        df_spK_Drop.loc[df_spK_Drop['spnam'].map(lambda x: x=='Methyl Mercury'),namK]=2.1
    cols=['spnam','CAS']+[namK]
    df_spK_Drop[cols].set_index('spnam').sort_values([namK]).to_csv('sp_'+namK+'.csv')
    print namK, len(df_spK_Drop)
    df_final=df_final.merge(df_spK_Drop,how='left', left_on='spnam', right_on='spnam')
    d.update({namK:Series(df_spK_Drop['spnam'].sort_values())})
cols=cols=['spnam','CAS']+nam_K
df_final[cols].set_index('spnam').sort_values('spnam').to_csv('sp_all.csv')
df_nam=DataFrame(d)
print df_nam

