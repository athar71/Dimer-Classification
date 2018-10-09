#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 15:32:03 2017

@author: athar
"""
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score
import numpy as np
#from math import floor


#==============================================================================
features=np.loadtxt('features_train_ssdu.txt')
label=np.loadtxt('label_train_ssdu.txt')
#==============================================================================
#==============================================================================
#features_test=np.loadtxt('features_test.txt')
#label_test=np.loadtxt('label_test.txt')
#features_train=np.loadtxt('features_train.txt')
#label_train=np.loadtxt('label_train.txt')
#==============================================================================
#label_test=-1*label_test
#==============================================================================
#features=np.concatenate((features_train, features_test), axis=0)
#label=np.concatenate((label_train, label_test), axis=0)
 
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.4, random_state=0)
train_index, test_index=sss.split(features,label).__next__()
X_test=features[test_index]
Y_test=label[test_index]
X_train=features[train_index]
Y_train=label[train_index]
# #==============================================================================
#==============================================================================
# X_test=np.loadtxt('features_test.txt')
# Y_test=np.loadtxt('label_test.txt')
# X_train=np.loadtxt('features_train.txt')
# Y_train=np.loadtxt('label_train.txt')
#==============================================================================

#create and train the random forest
#multi-core CPUs can use: rf = RandomForestClassifier(n_estimators=100, n_jobs=2)
rf = RandomForestClassifier(n_estimators=1000)
rf.fit(X_train, Y_train)

y_pred = rf.predict_proba(X_test)
#AUC_tree=roc_auc_score(Y_test==1, y_pred[:, 1])
AUC_tree=roc_auc_score(Y_test==1, y_pred[:, 1])
score_tree=rf.score(X_test,Y_test)

   