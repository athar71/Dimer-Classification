#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 10:10:11 2017

@author: athar
"""

from sklearn import svm
from sklearn.linear_model import LogisticRegression
import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedShuffleSplit
#X = [[0, 0], [1, 1]]
#y = [0, 1]
#==============================================================================
#==============================================================================
# X_test=np.loadtxt('features_test.txt')
# Y_test=np.loadtxt('label_test.txt')
# X_train=np.loadtxt('features_train.txt')
# Y_train=np.loadtxt('label_train.txt')
#  
#==============================================================================
#==============================================================================
#==============================================================================
#features_test=np.loadtxt('features_test.txt')
#label_test=np.loadtxt('label_test.txt')
 #label_test=-1*label_test
features_train=np.loadtxt('features_train_ssdu.txt')
label_train=np.loadtxt('label_train_ssdu.txt')
#features=np.concatenate((features_train, features_test), axis=0)
#label=np.concatenate((label_train, label_test), axis=0)
 
sss = StratifiedShuffleSplit(n_splits=1, test_size=0.4, random_state=0)
train_index, test_index=sss.split(features_train,label_train).__next__()
X_test=features_train[test_index]
Y_test=label_train[test_index]
X_train=features_train[train_index]
Y_train=label_train[train_index]
#==============================================================================
#clf = svm.SVC(probability=True)
# clf = svm.SVC(probability=True)
#clf = LogisticRegression()

clf=svm.SVC(probability=True, kernel='linear')
clf.fit(X_train, Y_train)  

#clf = svm.SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0,
#    decision_function_shape=None, degree=3, gamma='auto', kernel='rbf',
#    max_iter=-1, probability=False, random_state=None, shrinking=True,
#    tol=0.001, verbose=False)
#clf.fit(X_train, Y_train)  
#==============================================================================
#clf = svm.SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0, degree=3, kernel='rbf', max_iter=-1, probability=True, random_state=None,
# shrinking=True, tol=0.001, verbose=False)
#clf.fit(X_train, Y_train)  
#==============================================================================
y_pred = clf.predict_proba(X_test)
AUC_svm=roc_auc_score(Y_test==1, y_pred[:, 1])
score_svm=clf.score(X_test,Y_test)