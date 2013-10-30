import numpy as np
import pylab as pl
from sklearn import svm, datasets

# import some data to play with
iris = datasets.load_iris()
X = iris.data[:, :2]  # we only take the first two features. We could
                      # avoid this ugly slicing by using a two-dim dataset
Y = iris.target

x1 = iris.data[:, :1]

h = .02  # step size in the mesh

# we create an instance of SVM and fit out data. We do not scale our
# data since we want to plot the support vectors
C = 1.0  # SVM regularization parameter
svc = svm.SVC(kernel='linear', C=C).fit(X, Y)
rbf_svc = svm.SVC(kernel='rbf', gamma=0.7, C=C).fit(X, Y)
poly_svc = svm.SVC(kernel='poly', degree=3, C=C).fit(X, Y)
lin_svc = svm.LinearSVC(C=C).fit(X, Y)

# create a mesh to plot in
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                     np.arange(y_min, y_max, h))

# title for the plots
titles = ['SVC with linear kernel',
          'SVC with RBF kernel',
          'SVC with polynomial (degree 3) kernel',
          'LinearSVC (linear kernel)']


for i, clf in enumerate((svc, rbf_svc, poly_svc, lin_svc)):
    # Plot the decision boundary. For that, we will assign a color to each
    # point in the mesh [x_min, m_max]x[y_min, y_max].
    pl.subplot(2, 2, i + 1)
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])
    irisPred = clf.predict(np.c_[6, 3,2])
    irisPred = clf.predict(X)

    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    pl.contourf(xx, yy, Z, cmap=pl.cm.Paired)
    pl.axis('off')

    # Plot also the training points
    pl.scatter(X[:, 0], X[:, 1], c=Y, cmap=pl.cm.Paired)

    pl.title(titles[i])

pl.show()

  
### re-do example with 3 dim
X = iris.data[:, :3]  
Y = iris.target
# we create an instance of SVM and fit out data. We do not scale our
# data since we want to plot the support vectors
C = 1.0  # SVM regularization parameter
svc = svm.SVC(kernel='linear', C=C).fit(X, Y)
rbf_svc = svm.SVC(kernel='rbf', gamma=0.7, C=C).fit(X, Y)
poly_svc = svm.SVC(kernel='poly', degree=3, C=C).fit(X, Y)
lin_svc = svm.LinearSVC(C=C).fit(X, Y)

# create a mesh to plot in
x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
z_min, z_max = X[:, 2].min() - 1, X[:, 2].max() + 1
xx, yy, zz = np.meshgrid(np.arange(x_min, x_max, h),
                     np.arange(y_min, y_max, h),
                     np.arange(z_min, z_max, h))
# title for the plots
titles = ['SVC with linear kernel',
          'SVC with RBF kernel',
          'SVC with polynomial (degree 3) kernel',
          'LinearSVC (linear kernel)']
for i, clf in enumerate((svc, rbf_svc, poly_svc, lin_svc)):
    # Plot the decision boundary. For that, we will assign a color to each
    # point in the mesh [x_min, m_max]x[y_min, y_max].
    pl.subplot(2, 2, i + 1)
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel(),zz.ravel()])
    # irisPred = clf.predict(np.c_[6, 3,2]) #random set of values in the space
    irisPred = clf.predict(X)
    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    pl.contourf(xx[:,:,1], yy[:,:,1], Z[:,:,1], cmap=pl.cm.Paired)
    pl.axis('off')
    # Plot also the training points
    pl.scatter(X[:, 0], X[:, 1], c=Y, cmap=pl.cm.Paired)
    pl.title(titles[i])
pl.show()



### example 2  
# we create 40 separable points
np.random.seed(0)
X = np.r_[np.random.randn(20, 2) - [2, 2], np.random.randn(20, 2) + [2, 2]]
Y = [0] * 20 + [1] * 20

# fit the model
clf = svm.SVC(kernel='linear')
clf.fit(X, Y)

# get the separating hyperplane
w = clf.coef_[0]
a = -w[0] / w[1]
xx = np.linspace(-5, 5)
yy = a * xx - (clf.intercept_[0]) / w[1]

# plot the parallels to the separating hyperplane that pass through the
# support vectors
b = clf.support_vectors_[0]
yy_down = a * xx + (b[1] - a * b[0])
b = clf.support_vectors_[-1]
yy_up = a * xx + (b[1] - a * b[0])

# plot the line, the points, and the nearest vectors to the plane
pl.plot(xx, yy, 'k-')
pl.plot(xx, yy_down, 'k--')
pl.plot(xx, yy_up, 'k--')

pl.scatter(clf.support_vectors_[:, 0], clf.support_vectors_[:, 1],
           s=80, facecolors='none')
pl.scatter(X[:, 0], X[:, 1], c=Y, cmap=pl.cm.Paired)

pl.axis('tight')
pl.show()



## load in cmap data for svm
# two common drugs - 
CM = mu.CMapMongo()
goldQuery = CM.find({'is_gold' : True}, #, 
        {'sig_id':True,'pert_id':True,'cell_id':True,'pert_time':True,'is_gold':True,'pert_iname':True},
        toDataFrame=True)
grped = goldQuery.groupby(['pert_id','cell_id'])
grpSize = grped.size()
grpSize = grpSize.order(ascending=False)
grp1 = grpSize.index[0]
grp2 = grpSize.index[1]
igrp1 = grped.groups[grp1]
igrp2 = grped.groups[grp2]
grp1Quer = goldQuery.ix[igrp1]
grp2Quer = goldQuery.ix[igrp2]
grp1Sigs = grp1Quer['sig_id']
grp2Sigs = grp2Quer['sig_id']

sigList = list(grp1Sigs)
sigList.extend(grp2Sigs)
#load in expression data for the two sets of signatures
afPath = cmap.score_path
gt = gct.GCT()
gt.read(src=afPath,cid=sigList,rid='lm_epsilon')
zFrm = gt.frame
zFrm = zFrm.T
probeIDs = zFrm.columns
## asign labels
pert_ids = [x.split(':')[1][:13] for x in zFrm.index]
# create labels
labels = []
for brd in pert_ids:
  if brd == 'BRD-A19500257':
    labels.append(1)
  else:
    labels.append(0)
labels = np.array(labels)
zFrm['labels'] = labels

### perform leave one out validation
predictDict = {}
for sig in zFrm.index:
  droppedFrm = zFrm[zFrm.index != sig] # remove test signature from training
  trainFrm = droppedFrm.reindex(columns=probeIDs)
  labelsTrain = droppedFrm['labels'].values
  C = 1.0  # SVM regularization parameter
  svc = svm.SVC(kernel='linear', C=C).fit(trainFrm.values, labelsTrain)
  zTest = zFrm.ix[sig,probeIDs]
  linPred = svc.predict(zTest.values)
  predictDict[sig] = linPred[0]
predSer = pd.Series(predictDict)
predSer.name = 'svm_prediction'
zFrm = pd.concat([zFrm,pd.DataFrame(predSer)],axis=1)
realVPred = zFrm.ix[:,['labels','svm_prediction']]
accuracyArray = zFrm['labels'] == zFrm['svm_prediction']
accuracyRate = accuracyArray.sum()/float(accuracyArray.shape[0])

## make dataframe
# arrayList = [zFrm.values[ix,:] for ix in np.arange(zFrm.shape[0])]
# pd.DataFrame(arrayList,index=zFrm.index,columns=['z_scores'])
# train



