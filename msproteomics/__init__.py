# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash, html, dcc
from flask import Flask
import plotly.express as px
import pandas as pd
import numpy as np 
import os
PNUM = 'Ecoli TMT'

def rename_ratios4(dfIn):
    """
    Adapted for Peptide Groups table
    """
    renaming_dict = {}
    for c in dfIn.columns:
        if 'Abundance Ratio' in c:
            new_name = c.split(' ')[-3] + '/' + c.split(' ')[-1]
            renaming_dict[c] = new_name
        elif 'Abundances Grouped' in c:
            new_name = c.split(' ')[-1]
            renaming_dict[c] = new_name
        renaming_dict['Master Protein Accessions'] = 'Accession'
    dfOut = dfIn.rename(renaming_dict,axis='columns')
    return dfOut


def variability_groupby(df, ratio_columns, accession_col='Accession',
                          use_peptides='all',shared_keyword='Shared', unique_infocol='Quan Info'):
    if use_peptides == 'unique':
        df = df[~df[unique_infocol].isin((shared_keyword,))].copy()
    elif use_peptides == 'all':
        pass
    else:
        df = None
        
    dfN = df[ ([accession_col,] + ratio_columns) ]
    
    dfMeans = dfN.groupby([accession_col]).mean()
    dfDev = dfN.groupby([accession_col]).std()
    dfVar = np.divide(dfDev, dfMeans)
    dfVar = dfVar*100
    dfVar.replace(0, np.nan, inplace=True)
    dfVar = dfVar.round(2)
    
    return dfVar

df5 = pd.read_csv(
    os.path.join('data',
        "PXD007647_Reproc_TMT-set-2_8fracs_PeptideGroups.txt"),
    sep='\t'
)



df5 = rename_ratios4(df5)


#server = Flask(__name__)
#server.config['SECRET_KEY'] = 'AdrianAlicia20222021249k4hvd'
app = Dash(__name__)

dfPeptVar = variability_groupby(
    df5[ (['Accession', ] + list(df5.columns[15:24]) ) ].dropna(axis='rows'),
    ratio_columns = list( df5.columns[15:24] )
)


fig1 = px.box( pd.DataFrame(np.log10(
    dfPeptVar.replace( 0, np.nan))),
             labels={
                     "value": "Log10 Peptide Variability (%)",
                     #"variable" : "abcd",
                 },
    title=f'{PNUM} Peptide Variability (%) Within Proteins'
   )



#fig2 = px.density_heatmap(data_frame= pd.DataFrame(np.log2(df5.iloc[:, 15:24].dropna(axis='rows'))).corr(), 
fig2 = px.imshow(pd.DataFrame(np.log2(df5.iloc[:, 15:24].dropna(axis='rows')).corr().round(2)),
                text_auto = True,
                color_continuous_scale='RdBu_r',
                title = f'{PNUM} Pearson Correlations on Peptide Level'
                )



fig3 =  px.scatter_matrix(data_frame=pd.DataFrame(np.log2(df5.iloc[:, 15:24].dropna(axis='rows'))),
                              title= f'{PNUM} Pair Plots on Peptide Level' , 
                              
                               )
# fig3.update_xaxes(tickangle= 45)

""" 
df5['Cys_Peptide'] = [
    'Yes' if '[C' in x else 'No' for x in df5['Modifications']
]
fig4 = dashbio.Clustergram(pd.DataFrame(np.log2(
        df5[
            df5['Cys_Peptide'] == 'No'
        ].iloc[
            :, 15:24
        ].dropna(axis='rows')
    )))
#, title = f'{PNUM} Modified Cys Peptide Log2 Relative Ints'
"""

app.layout = html.Div(children=[
    html.H1(children='MS Proteomic Data Analysis', style={'text-align': 'center','color': 'lightgray', 'padding-top': '70px',
  'padding-bottom': '70px','background-color': 'black', 'margin-top':'0px', 'margin-bottom':'0px'}),

    html.Div(children= html.H2(f'{PNUM} Peptide Variability (%) Within Proteins',style={'text-align': 'center' , 'background-color': 'lightgrey',
  'color': 'blue', 'padding-top': '50px',
  'padding-bottom': '50px','margin-top':'0px'})),
    html.Div(children=[
    dcc.Graph(
        id='graph1',
        figure=fig1
    ),
    html.Div(children= html.H2(f'{PNUM} Pearson Correlations on Peptide Level', style={'text-align': 'center' , 'background-color': 'lightgrey',
  'color': 'blue', 'padding-top': '50px',
  'padding-bottom': '50px',}))]),

    dcc.Graph(
        id='graph2',
        figure=fig2
    ),
    html.Div(children=[
        html.Div(children=html.H2(f'{PNUM} Pair Plots on Peptide Level',
    ),style={'text-align': 'center' , 'background-color': 'lightgrey',
  'color': 'blue', 'padding-top': '50px',
  'padding-bottom': '50px',}
    ),

    dcc.Graph(
        id='graph3',
        figure=fig3, style={'width': '90vw', 'height': '90vh','text-align': 'center', 'margin': 'auto'}
    )
    ],)
    ,
      
    
    ])


