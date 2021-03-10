#function ARGUMENTs
nodalCoordinatesThisRun = [0,0,0,0,0,3]
axialLoadThisRun=[0,0,-100,0,0,0] #apply -100kN in the direction of gravity



from singleColumnFunctionFile import singleColumnAnalysisFunction

data=singleColumnAnalysisFunction(nodalCoordinatesThisRun,axialLoadThisRun) #3 element arrat containing the returned variables (node reactions, node displacement and element forces )

print("done")