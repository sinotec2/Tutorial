# coding=utf-8
from pandas import *
import numpy as np

def echos():
    print 'input xls total records for sheet1:',len(df1)
    print 'input xls total records for sheet2:',len(df2)
    print 'input xls total records for ',wsname,':',len(df_cat)
    print 'sheet1 cols:',df1.columns
    print 'sheet2 cols:'
    for i in df2.columns:print i
    print wsname,' cols:',df_cat.columns
    plt=list(df1[u'管制編號'])
    plant=set(plt)
    print 'num. of plants selected:',len(plant_s)
    print 'out of all plants:',len(plant)
    print 'the water pollution sources in taiwan:'
    #print pivot_table(df_all,index=["cat_all"],values=["reg_all"],aggfunc='count')    
    print 'the water pollution sources selected:'
    #print pivot_table(df_s,index=["category_selected"],values=["plant"],aggfunc='count')    
    print 'taiwan 104 yr COD sources and total COD(Ton)'
    print len(list(df_2pv['COD'])),sum(list(df_2pv['COD']))
    print 'merging the sheet1_selected and sheet2'
    print df_pv
    print 'categorized COD sum checking'
    print 'sum check', sum(list(df_pv['COD']))
    return

#read the input file
#fname=u'160617  國家溫室氣體清冊需求資料-水保處提供.xlsx'
fname=u'160617.xlsx'
wsname=u'11業'
df1 = read_excel(fname,sheetname='sheet1',skiprows=1, parse_cols = 'A,G,J')
df2 = read_excel(fname,sheetname='sheet2',skiprows=1, parse_cols = 'A:G',na_values=['NA'], convert_float=True)
df_cat = read_excel(fname,sheetname=wsname,skiprows=0, parse_cols = 'A:B')
print 'the input file is read'
#first of all, reduce the whole database and accesslorate the process
(cat,plt,reg)=([],[],[])
#using set() but not list() to avoid duplicates(different facilities)
for p in set(list(df1[u'管制編號'])):
    boo=df1[u'管制編號']==p   #the boolean criteria
    cat.append(list(df1[boo][u'行業別'])[0][6:]) #store the cat. without the number
    reg.append(list(df1[boo][u'區域類別(如是否位於工業區)'])[0])
    plt.append(p)
a={'plant_all':Series(plt),'cat_all':Series(cat),'reg_all':Series(reg)}
df_all=DataFrame(a).sort_values('plant_all')    #the new, reduced df of all sources
#secondly, do the section
cat_nam=list(df_cat[u'業別'])
reg_nam=list(df_cat[u'工業區'])
(categ_s,plant_s)=([],[])
for c in cat_nam:
    for r in reg_nam[0:2]: #other are blank (NaN)
        #boolean for matching categories and regions
        #categories with number must be skipped(already cut by previous step)
        boo=(df_all['cat_all']== c)&(df_all['reg_all']==r)    #multipled criteria
        l=list(df_all[boo]['cat_all'])
        if len(l)>=1:    #maybe no item at all
            for p in set(list(df_all[boo]['plant_all'])): 
                boo=df_all['plant_all']==p
                #store the first record is enough and only
                categ_s.append(list(df_all[boo]['cat_all'])[0]) 
                plant_s.append(p)        
a={'plant':Series(plant_s),'category_selected':Series(categ_s)}
df_s=DataFrame(a).sort_values('plant') #sources with certain categories and regions
print 'sheet1 is done, now dealwith sheet2, the COD'
#reading the COD for each plant and season
#screening the year and selecting the plants
boo=(df2[u'管制編號'].map(lambda x: x in plant_s ))
df2_s=df2[boo]
boo=(df2_s[u'申報區間(起)'].map(lambda x: '104' in x)) 
df2_yr=df2_s[boo]
#calculate the COD emissions                                                     
COD=((df2_yr[u'進流污水量(噸)'] * df2_yr[u'進流口COD值(mg/L)'] - \
   df2_yr[u'放流水量(噸)'] * df2_yr[u'放流口COD值(mg/L)'] ) /10**6).clip(0)
a={'plant':Series(df2_yr[u'管制編號']),'COD':Series(COD)}
df_2=DataFrame(a).sort_values('plant')
df_2pv=pivot_table(df_2,index=["plant"],values=["COD"],aggfunc=np.sum)
df_2pv['plant']=df_2pv.index #reset the plant as a column
print 'COD has been done, now doing the merging'
a=set(plant_s)-set(df_2pv['plant']) #drop the plant without COD data
for i in a:
    df_s.drop(df_s[df_s['plant']==i].index, inplace=True)
df=merge(df_s,df_2pv,on=['plant'])
df_pv=pivot_table(df,index=["category_selected"],values=["COD"],aggfunc=np.sum)
writer = ExcelWriter('output.xlsx')
df_pv.to_excel(writer,'Sheet1')
df_s.to_excel(writer,'Sheet2')
df_2pv.to_excel(writer,'Sheet3')
writer.save()

echos()
