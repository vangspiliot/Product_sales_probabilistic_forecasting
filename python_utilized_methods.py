from sklearn import ensemble
import numpy as np

def X_Matrix(insample, window):              #form the X matrix for the training of the ml model
    
    n = len(insample)
    if (n > window):
    
        d = n - window
        tr = np.empty((d, window))
    
        for i in range(d):
            for j in range(window):
                tr[i][j] = insample[i+j]
    else:
        value = insample[0]
        new_vect = []
        for i in range(window-len(insample)+1):
            new_vect.append(value)
        for i in range(len(insample)-1):
            new_vect.append(insample[i])
        tr = np.array(new_vect)
        tr = np.reshape(tr, (1,window))
        

    return tr


def y_Matrix(insample, window):              #form the y matrix for the training of the ml model
    
    n = len(insample)
    if (n > window):

        d = n - window
        tr = []

        for i in range(d):
            tr.append(insample[i+window])
    else:
        tr = insample[-1]

    F = np.array(tr)

    return F


def setforprediction(insample, window):     #form the x matrix for predicting using the ml model
    
    wanted = np.array(insample[-window:])
    wanted = np.reshape(wanted, (1,window))
    
    return wanted


def gb_calc(insample, horizon, quant):  #predict quantile using a gradient boosting model
    
    #insample: insample of a given time series
    #horizon: forecasting horizon
    #quant: wanted quantile

    window = 28                                         #window for forming the matrices for ml model training
    estimators = 500                                    #No of estimators for ml model
    lr = 0.01                                           #learning rate for the gradient boosting model
    ss = 10                                             #the minimum number of samples required to split an internal node

    x_train = X_Matrix(insample, window)
    y_train = y_Matrix(insample, window)
    x_test = setforprediction(insample, window)

    gbf = ensemble.GradientBoostingRegressor(   loss='quantile', alpha=quant, n_estimators=estimators, 
                                                max_depth=3, learning_rate=lr, min_samples_split=ss)
    
    gbf.fit(x_train, y_train)                           #train

    predicted_quantile = gbf.predict(x_test)            #predict

    if float(predicted_quantile) <=0:
        predicted_quantile = 0
    
    return float(predicted_quantile)


def rf_calc(insample, horizon, quant):  #predict quantile using a random forest model
    
    #insample: insample of a given time series
    #horizon: forecasting horizon
    #quant: wanted quantile

    window = 28                                             #window for forming the matrices for ml model training
    estimators = 500                                        #No of estimators for ml model
    ss = 10                                                 #the minimum number of samples required to split an internal node

    x_train = X_Matrix(insample, window)
    y_train = y_Matrix(insample, window)
    x_test = setforprediction(insample, window)

    
    rf = ensemble.RandomForestRegressor(n_estimators=estimators, 
                                min_samples_leaf=2, min_samples_split=ss, random_state=3, 
                                verbose=True, 
                                n_jobs=-1)                  # Use maximum number of cores.

    rf.fit(x_train, y_train)                                #train

    rf_preds = []
    
    for estimator in rf.estimators_:
        rf_preds.append(estimator.predict(x_test))

    rf_preds = np.array(rf_preds).transpose()

    predicted_quantile = np.percentile(rf_preds, quant * 100, axis=1)

    if float(predicted_quantile) <=0:
        predicted_quantile = 0
    
    return float(predicted_quantile)


import statsmodels.api as sm


def lqr_calc(insample, horizon, quant): #predict quantile using a linear quantile regression model

    #insample: insample of a given time series
    #horizon: forecasting horizon
    #quant: wanted quantile

    Y_train = insample                                                              #create the lr model where, Y_train is the observations of the insample and X_train a vestor representing time
    X_train = range(1,(len(insample)+1))
    X_train = sm.add_constant(X_train)
    model = sm.QuantReg(Y_train,X_train)

    X_pred = range((len(insample)+1),(len(insample)+horizon+1))                     #set the wanted time vector for the forecasting horizon
    X_pred = sm.add_constant(X_pred)

    predicted_quantiles = model.fit(q=quant).predict(X_pred)                        #predict

    for i in range(len(predicted_quantiles)):
        if predicted_quantiles[i] <=0:
            predicted_quantiles[i] = 0
    
    return predicted_quantiles

