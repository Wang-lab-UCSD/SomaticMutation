import os, sys
import numpy as np
import math
from scipy import stats
import pandas
from keras.utils import plot_model
from keras.models import Sequential
from keras.models import Model
from keras.layers import Input
from keras.layers import Dense
from keras.layers import Dropout
from keras.constraints import maxnorm
from keras import regularizers
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import matplotlib.pyplot as plt
import tensorflow as tf
plt.switch_backend('agg')
from keras.constraints import nonneg
from keras import backend as K
from keras.models import load_model
from keras.utils import plot_model
import keras
import gc 

def correlation_coefficient_loss(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = K.mean(x)
    my = K.mean(y)
    xm, ym = x-mx, y-my
    r_num = K.sum(tf.multiply(xm,ym))
    r_den = K.sqrt(tf.multiply(K.sum(K.square(xm)), K.sum(K.square(ym))))
    r = r_num / r_den
    r = K.maximum(K.minimum(r, 1.0), 0)
    return 1 - K.square(r)

def likelihood_loss(y_true,y_pred):
    x = y_true/y_pred
    y = K.log(y_pred+1)
    C_BS = K.sum(x+y)
    C_BS = 0.5*C_BS
    return C_BS

# huber loss
def huber(true, pred, delta=1):
    tmp = np.where(np.abs(true-pred) < delta , 0.5*((true-pred)**2), delta*np.abs(true - pred) - 0.5*(delta**2))
    return np.sum(tmp)

def huber_loss(y_true, y_pred):
    return tf.losses.huber_loss(y_true,y_pred,delta=0.5)

def logcosh(true, pred):
    loss = np.log(np.cosh(pred - true))
    return np.mean(loss)

def correlation_coefficient(y_true, y_pred):
    x = y_true
    y = y_pred
    mx = K.mean(x)
    my = K.mean(y)
    xm, ym = x-mx, y-my
    r_num = K.sum(tf.multiply(xm,ym))
    r_den = K.sqrt(tf.multiply(K.sum(K.square(xm)), K.sum(K.square(ym))))
    r = r_num / r_den
    r = K.maximum(K.minimum(r, 1.0), -1.0)
    return r

def error_model():
    model = Sequential()
    model.add(Dense(Num_Feature, input_dim=Num_Feature, kernel_initializer='random_uniform', activation='relu',kernel_constraint=maxnorm(2)))
    model.add(Dropout(0.2))
    model.add(Dense(Num_Feature, kernel_initializer='random_uniform', activation='relu',kernel_constraint=maxnorm(5)))
    model.add(Dropout(0.1))
    model.add(Dense(1, kernel_initializer='random_uniform',activation='exponential'))
    model.compile(loss=likelihood_loss, optimizer='adam',metrics=[correlation_coefficient])
    return model


def tilted_loss(y,f):
    q = 0.75
    e = (y-f)
    return K.mean(K.maximum(q*e, (q-1)*e), axis=-1)

def load_all_models(n_models_Seq,Dir,Name):
    all_models = list()
    for i in n_models_Seq:
        # define filename for this ensemble
        filename = Dir+ Name +str(i+1) + '.h5'
        # load model from file
        #model = KerasRegressor(build_fn=baseline_model, epochs=500, batch_size=2000, verbose=1)
        model = load_model(filename, custom_objects={'tilted_loss': tilted_loss,'correlation_coefficient': correlation_coefficient})
        # add to list of members
        all_models.append(model)
        print('>loaded %s' % filename)
    return all_models


def load_contribution_models(n_models_Seq,Dir):
    all_models = list()
    for i in n_models_Seq:
        # define filename for this ensemble
        filename = Dir+ 'model_Contribution' +str(i+1) + '.h5'
        # load model from file
        #model = KerasRegressor(build_fn=baseline_model, epochs=500, batch_size=2000, verbose=1)
        model = load_model(filename, custom_objects={'correlation_coefficient': correlation_coefficient})
        # add to list of members
        all_models.append(model)
        print('>loaded %s' % filename)
    return all_models

def load_Var_models(n_models,Dir):
    all_models = list()
    for i in range(n_models):
        # define filename for this ensemble
        filename = Dir+ 'model_' +str(i + 1) + '.h5'
        # load model from file
        model = KerasRegressor(build_fn=baseline_model, epochs=500, batch_size=2000, verbose=1)
        model.model = load_model(filename, custom_objects={'likelihood_loss': likelihood_loss,'correlation_coefficient': correlation_coefficient})
        # add to list of members
        all_models.append(model)
        print('>loaded %s' % filename)
    return all_models



def stacked_dataset(members, inputX):
    N = len(members)
    stackX = np.zeros( (len(inputX),N) )
    for i in range(N):
        model = members[i]
        yhat = model.predict(inputX, verbose=0)
        yhat = yhat.reshape((len(inputX),))
        stackX[:,i] = yhat
    return stackX

def stacked_dataset_multipleInput(members, inputX):
    N = len(members)
    stackX = np.zeros( (len(inputX[0]),N) )
    for i in range(N):
        model = members[i]
        yhat = model.predict(inputX, verbose=0)
        yhat = yhat.reshape((len(inputX[0]),))
        stackX[:,i] = yhat
    return stackX



def ContributionMean_multipleInput(members, inputX):
    N = len(members)
    stackX = np.zeros((len(inputX),inputX.shape[1]+1))
    Test_X1 = np.ones((len(inputX),1))
    Test_X1 = np.append(inputX, Test_X1, axis=1)
    for i in range(N):
        model = members[i]
        yhat = model.predict(inputX, verbose=0)
        yhat = Test_X1*yhat
        stackX = stackX+yhat
    stackX = stackX/N
    return stackX


def ContributionMean(members, inputX):
    N = len(members)
    stackX = np.zeros((len(inputX),inputX.shape[1]))
    for i in range(N):
        model = members[i]
        yhat = model.predict(inputX, verbose=0)
        yhat = inputX*yhat
        stackX = stackX+yhat
    stackX = stackX/N
    return stackX

def ContributionWeight(members, inputX):
    N = len(members)
    stackX = np.zeros((len(inputX),inputX.shape[1]))
    for i in range(N):
        model = members[i]
        yhat = model.predict(inputX, verbose=0)
        stackX = stackX+yhat
    stackX = stackX/N
    return stackX





def baseline_model(Num_Feature):
 visible1 = Input(shape=(Num_Feature,),name="Input1")
 hidden1 = Dense(int(Num_Feature/2), kernel_initializer='random_uniform', activation='relu',kernel_constraint=maxnorm(2))(visible1)
 hidden2 = Dropout(0.01)(hidden1)
 hidden3 = Dense(int(Num_Feature/10-20), kernel_initializer='random_uniform', activation='relu',kernel_constraint=maxnorm(5))(hidden2)
 hidden4 = Dropout(0.01)(hidden3)
 hidden5 = Dense(int(Num_Feature/5-30), kernel_initializer='random_uniform', activation='relu',kernel_constraint=maxnorm(5))(hidden4)
 hidden6 = Dropout(0.1)(hidden5)
 Contribution = Dense(Num_Feature, kernel_initializer='random_uniform',activation='linear',name='Contribution')(hidden6)
 output = keras.layers.dot([visible1,Contribution], axes=1,normalize=False)
 #Bias = BIAS(1,kernel_initializer='random_uniform', activation='relu')
 ContextualRegression = Model(inputs=visible1, outputs=output)
 adam = keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=None, decay=0.001, amsgrad=False)
 ContextualRegression.compile(loss='mean_squared_error', optimizer=adam,metrics=[correlation_coefficient])
 print(ContextualRegression.summary())
 return ContextualRegression



#----Start here--------------------------
DataFile = sys.argv[1]
Out = sys.argv[2]
GPU = sys.argv[3]
Epoch = int(sys.argv[4])
ScriptFolder = sys.argv[5]

Out=Out+'Normalize0-1_cluterMotif_Zscore_noDHS'

os.environ["CUDA_VISIBLE_DEVICES"] = str(GPU)
config = tf.ConfigProto()
config.gpu_options.allow_growth = True  # dynamically grow the memory used on the GPU
config.log_device_placement = True  # to log device placement (on which device the operation ran)
sess = tf.Session(config=config)


NodeSelf = pandas.read_csv(DataFile,sep="\t")

Zero_ID = NodeSelf[NodeSelf.columns[3:len(NodeSelf.columns)]]>0
COL_SUM = Zero_ID.apply(lambda x: x.sum(), axis=1)

NodeSelf = NodeSelf[COL_SUM>3]
NodeSelf.index = range(0,len(NodeSelf))


TotalNum = NodeSelf['Num_mutation'].sum()/1000000.0
Y = NodeSelf['Num_mutation']/NodeSelf['L']*1000.0/TotalNum
Y = np.log2(Y+1)
Scale = Y.mean()
STD = np.sqrt(Y.var())
Y = (Y-Scale)/STD


f = plt.figure()
plt.hist(Y, bins = 200)
plt.xlabel('Z-score(log2(MutationRate+1))')
plt.show()
f.savefig(Out+"/HistgramOfZscoreOfMutationRate.pdf", bbox_inches='tight')
plt.close()




Num = len(NodeSelf.columns)

Num_Feature = len(NodeSelf.columns)-3


Train = NodeSelf[NodeSelf.columns[3:Num]].values

Feature = list(NodeSelf.columns[3:Num])


del NodeSelf
gc.collect()


#Train = np.log2(Train+1)

MAX_col = np.max(np.abs(Train),axis=0)
tmp = np.where(MAX_col==0)[0]
if len(tmp)>0:
    MAX_col[tmp] = 1



Train /=  MAX_col


B = 10
L = range(len(Train))
kf = KFold(n_splits=B,shuffle=True)
kf.get_n_splits(L)




for h in ['0','0.1']:
    if h=='0':
        tmp = Feature
    else:
        Cluster = pandas.read_csv(ScriptFolder+'/Cluster2431Motif_Info_h'+h+'.txt',sep="\t")
        tmp = list(Cluster['ID'])+['DHS']
    Feature_Index = np.where(np.isin(Feature,tmp))[0].tolist()
    count = 0
    Result = pandas.DataFrame(np.zeros(2*B).reshape((B,2)),index=range(0,B),columns=['CR_NodeSelf_train','CR_NodeSelf_test'])
    for train_index, test_index in kf.split(L):
        #---NodeSelf--------------------
        Train_X = Train[train_index,:]
        Train_Y = Y[train_index]
        Test_X = Train[test_index,:]
        Test_Y = Y[test_index]
        Train_X = Train_X[:,Feature_Index]
        Test_X = Test_X[:,Feature_Index]
        ####---CR NodeSelf train----------------------
        ContextualRegression = baseline_model(Train_X.shape[1])
        ContextualRegression.fit(Train_X, Train_Y, epochs=Epoch, batch_size=500, verbose=True)
        filename = Out+'/CR_model_NodeSelf_MotifDHS_Epoch'+str(Epoch) + '_h'+h+'_'+str(count + 1) + '.h5'
        ContextualRegression.save(filename)
        layer_name = 'Contribution'
        Contribuiton_layer_model = Model(inputs=ContextualRegression.input,outputs=ContextualRegression.get_layer(layer_name).output)
        Contribuiton_layer_model.compile
        filename = Out+'/CR_model_NodeSelf_MotifDHS_ContributionModel_Epoch'+str(Epoch) + '_h'+h+'_'+str(count + 1) + '.h5'
        Contribuiton_layer_model.save(filename)
        prediction = ContextualRegression.predict(Train_X)
        prediction = prediction.reshape((len(Train_X),))
        cor = np.corrcoef(prediction,Train_Y)[0,1]
        cor = round(cor,5)
        Result['CR_NodeSelf_train'][count] = cor
        spr = round(stats.spearmanr(prediction, Train_Y)[0],5)
        MAE = round(abs(Train_Y-prediction).mean(),5)
        MSE = round((pow(Train_Y-prediction,2)).mean(),5)
        MAX = max(max(prediction),max(Train_Y))+0.2
        MIN = min(min(prediction),min(Train_Y))-0.2
        f = plt.figure()
        plt.scatter(x=prediction,y=Train_Y,s=2)
        plt.xlim(MIN,MAX)
        plt.ylim(MIN,MAX)
        plt.xlabel('prediction')
        plt.ylabel('Z-score(log2(MutationRate+1))')
        plt.plot([MIN, MAX], [MIN, MAX], ls="--", c=".3")
        plt.title('cor='+str(cor)+' MAE='+str(MAE)+' MSE='+str(MSE)+' spearman='+str(spr))
        plt.show()
        f.savefig(Out+"/ScatterPLotForCR_NodeSelf_train_Epoch"+str(Epoch) + '_h'+h+'_'+str(count+1)+".pdf", bbox_inches='tight')
        print('>Saved %s' % filename)
        plt.close()  
        #--CR Node Self test-------------
        prediction = ContextualRegression.predict(Test_X)
        prediction = prediction.reshape((len(Test_X),))
        cor = np.corrcoef(prediction,Test_Y)[0,1]
        cor = round(cor,5)
        Result['CR_NodeSelf_test'][count] = cor
        spr = round(stats.spearmanr(prediction, Test_Y)[0],5)
        MAE = round(abs(Test_Y-prediction).mean(),5)
        MSE = round((pow(Test_Y-prediction,2)).mean(),5)
        MAX = max(max(prediction),max(Test_Y))+0.2
        MIN = min(min(prediction),min(Test_Y))-0.2
        f = plt.figure()
        plt.scatter(x=prediction,y=Test_Y,s=2)
        plt.xlim(MIN,MAX)
        plt.ylim(MIN,MAX)
        plt.xlabel('prediction')
        plt.ylabel('Z-score(log2(MutationRate+1))')
        plt.plot([MIN, MAX], [MIN, MAX], ls="--", c=".3")
        plt.title('cor='+str(cor)+' MAE='+str(MAE)+' MSE='+str(MSE)+' spearman='+str(spr))
        plt.show()
        f.savefig(Out+"/ScatterPLotForCR_NodeSelf_test_Epoch"+str(Epoch) + '_h'+h+'_'+str(count+1)+".pdf", bbox_inches='tight')
        print('>Saved %s' % filename)    
        plt.close()
        count = count+1
        del Train_X,Test_X
        gc.collect()
    Result.to_csv(Out+'/predictionCor_Summary_CV'+str(B)+'Epoch'+str(Epoch) + '_h'+h+'.txt',index=0,sep='\t',header=True)   







