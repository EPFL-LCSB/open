import pandas as pd
import glob
from open.utils.postprocess import remove_log_calc_overallgamma
import time
'''
script is to merge h5 files generated for each sampled points
the combined data frames for the presented data in the manuscript can be found in /data folder 
in the repository
'''
def take_required_columns(df):
    for col in df.columns:
        #to remove binary variables from the main data
        #or col.startswith('z') \
        if col.startswith('d') or col.startswith('k_') \
                or col.startswith('m_') \
                or col.startswith('a')\
                or col.startswith('mt_') or col.startswith('t_')\
                or (('k' in col) and ('_' in col)):
            #split = str(col).split('_')
            #df[split[1]] = np.exp(df[col])
            df = df.drop(col, axis=1)
#todo make it general here it is assumed that if alpha there it is random bi-uni
   # if ('alpha' in df.columns):

    #    df['overall_gamma'] = (df.gamma_1 * df.gamma_3 * df.gamma_5 * df.gamma_6).round(3)
#and if not alpha it is uniuni 3 step so multiply three gammas -- not necessary just for post processing !
    #else:
     #   df['overall_gamma'] = (df.gamma_1 * df.gamma_2 * df.gamma_3 ).round(3)

    return df

#this is for the milp version postprocessing just remove the binary variables etc


def merge_files_milp(allFiles,percentage=1,alpha=0.5):
    list_ = []
    frame = pd.DataFrame()
    for file_ in allFiles:
        df = pd.read_hdf(file_)
        #a = file_.split('close_to_optimal_cutoff_')[1][0:3]
        #df['cut_off']=float(a)*np.ones(df.shape[0])
        #df=take_required_columns(df)
        #df=df[df.v>=(df.v.max())*percentage]
        list_.append(df)
    frame = pd.concat(list_, ignore_index=True)
    return frame



t=time.time()
'''provide path to your h5 files input'''
input_file='/open/projects/UniUni_wMILP/output_3step_fwd_170321_q_2.0/'
path=(input_file)
allFiles = glob.glob(path + "/*.h5")
n=100
chunks=[allFiles[i:i + n] for i in range(0, len(allFiles), n)]
list_=[]
for i in range(len(chunks)):
    allFiles=chunks[i]
    print(i)
#allFiles = glob.glob(path + "/*.h5")
    df=merge_files_milp(allFiles)
    #df['overall_gamma']=df['gamma_ov']
    list_.append(df)

frame = pd.concat(list_, ignore_index=True)

frame.to_hdf(input_file+'/df_combined.h5', key='s')

elapsed = time.time() - t
print('time for optim', elapsed)