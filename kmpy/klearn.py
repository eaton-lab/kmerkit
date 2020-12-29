#!/usr/bin/env python 

"""
A collection of scikit-learn machine learning applications for 
testing associations between genotype matrices (kmers)
and phenotypes using supervised learning, and for exploring 
patterns in kmer sharing among samples using unsupervised learning.


https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.RandomForestClassifier.html
https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.GradientBoostingClassifier.html#sklearn.ensemble.GradientBoostingClassifier

"""


import os
import sys
import numpy as np

import toyplot, toyplot.browser
from loguru import logger
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier



class Klearn:
    """
    ...
    """
    def __init__(self, name, workdir, phenos, trait):

        # store inputs
        self.name = name
        self.phenos = phenos
        self.trait = trait
        self.workdir = os.path.realpath(os.path.expanduser(workdir))
        self.prefix = os.path.join(self.workdir + f"_{self.name}")



    def load_pheno(self):
        """
        Load the phenotypes data and subsample taxa...
        """


    def load_matrix(self):
        """

        """
        # path to the genos matrix
        self.matrix_path = os.path.join(
            self.workdir + f"_kmatrix_{self.name}_var_genos.npy")

        # load matrix from disk as read-only
        self.matrix = np.memmap(
            filename=self.matrix_path,
            mode='r',
            dtype=np.bool,
        )
        self.matrix = self.matrix.reshape()


    def random_forest(self):
        """
        Random forest model to infer weights on combinations of kmers 
        that can best explain the observed phenotypes.
        """
        # TODO



    def rf_gradient_boost(self):
        """

        """

        # setup classifier
        params = {
            'n_estimators': 1200, 
            'max_depth': 3, 
            'subsample': 0.5,
            'learning_rate': 0.01, 
            'min_samples_leaf': 1, 
            'random_state': 3,
        }
        clf = GradientBoostingClassifier(**params)

        # split data into test/train
        # ...

        # fit training data
        clf.fit(X_train, y_train)
        
        # calculate score on test data
        acc = clf.score(X_test, y_test)




    def pca_samples(self, **kwargs):
        """
        Map samples into high dimensional space based on their 
        kmers. Runs a PCA first to reduce dimensionality.
        """
        mod_pca = PCA(**kwargs)
        res_pca = mod_pca.fit_transform(self.matrix[:, :500000])
        logger.info(f"PCA vars: {mod_pca.explained_variance_ratio_}")

        # plot result
        canvas = toyplot.Canvas(width=400, height=350)
        axes = canvas.cartesian()
        axes.scatterplot(
            res_pca[:, 0], 
            res_pca[:, 1],
            size=12,
            opacity=0.7,
            color=(
                toyplot.color.CategoricalMap().colors(self.phenos[self.trait])
            ),
        )
        axes.x.ticks.show = True
        axes.y.ticks.show = True
        toyplot.browser.show(canvas)
        return res_pca


    def tsne_samples(self, **kwargs):
        """
        t-distrubuted stochastic neighbor embedding
        """
        res = self.pca_samples()
        mod = TSNE(**kwargs)
        res2 = mod.fit_transform(res)

        # plot results
        canvas, _, _ = toyplot.scatterplot(
            res2[:, 0], 
            res2[:, 1],
            width=400,
            height=350,
            title=self.phenos.index,
            size=15, 
            opacity=0.7,
            color=(
                toyplot.color.CategoricalMap().colors(self.phenos[self.trait]),
            ),
        )
        toyplot.browser.show(canvas)



    def pca_kmers(self):
        """
        Map kmers into high dimensional space based on their 
        presence in samples.
        """
        #TODO

