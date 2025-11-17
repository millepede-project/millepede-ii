

def getParser():
    from argparse import ArgumentParser 
    ap = ArgumentParser("compareResults")
    ap.add_argument("--mode",type=str,choices=["res","log","mon"],default=None) 
    ap.add_argument("ref",type=str)
    ap.add_argument("test",type=str) 
    ap.add_argument("--lines","-l",type=int,default=-1) 
    return ap 

def parseResults(theFile, elements=-1):
    result = {}
    with open(theFile,"r") as fin: 
        for line in fin.readlines():
            tokens=line.split()
            try: 
                par = int(tokens[0])
                val = float(tokens[1])
                entries = int(tokens[4])
            except ValueError:  # skip documentation lines
                continue  
            except IndexError: # malformed lines - happens for underconstrained par
                 continue
            result[par]=(val,entries)
    return result
    
def parseResiduals(theFile, elements=-1):
    result = {}
    with open(theFile,"r") as fin: 
        for line in fin.readlines():
            tokens=line.split()
            try: 
                LFC = int(tokens[0])
                par = int(tokens[1])
                entries = int(tokens[2])
                median = float(tokens[3])
                rms = float(tokens[4])
                error = float(tokens[5])
            except ValueError:  # skip documentation lines
                continue  
            if not LFC in result.keys():
                result[LFC]={}
            result[LFC][par]=(entries,median,rms,error)
    return result

def parseLogFile(theFile):
    result = {}
    keys_iterTable = ["it" ,"fc","fcn_value", "dfcn_exp", " slpr", " costh" ] 
    keys_chi2 = ["Sum" ,"Chi^2","Ndf", "=" ] 
    with open(theFile,"r") as fin: 
        lines = fin.readlines()

        # Extract the lines describing the fit convergence 
        headerIndex = [i for i,l in enumerate(lines) if len([k for k in keys_iterTable if k in l])==len(keys_iterTable)]
        if len(headerIndex) == 1:
            convergence = []
            for thisLine in lines[headerIndex[0]+2:]:
                tokens = thisLine.split()
                if len(tokens)==0:
                    break # empty line at end of iteration printout 
                elif len(tokens) == 7 and int(tokens[0]) == 0:
                    # first iteration with inversion 
                    convergence += [(float(tokens[2]), 1.0, int(tokens[4]))]
                elif len(tokens) == 10 and int(tokens[0]) > 0:
                    # later iterations 
                    convergence += [(float(tokens[2]), float(tokens[3]), int(tokens[7]))] 
        c2Index = [i for i,l in enumerate(lines) if len([k for k in keys_chi2 if k in l])==len(keys_chi2)]
        if len(c2Index) == 1:
            for thisLine in lines[c2Index[0]:]:
                tokens = thisLine.split()
                if tokens[0] == "=":
                     result["Chi2"] = float(tokens[1])
                     break 
    result["convergence"] = convergence
    return result     

def compareResults(refFile,testFile,lines):
    refResult = parseResults(refFile,lines)
    testResult = parseResults(testFile,lines)
    success = True 

    # Step 1: Check if all labels appear in both files 
    uniqueToTest = [label for label in testResult.keys() if label not in refResult.keys()]
    uniqueToRef = [label for label in refResult.keys() if label not in testResult.keys()]
    if len(uniqueToTest) > 0:
            print("=== Difference: Some labels only appear in test file:")
            from pprint import pprint 
            pprint(uniqueToTest) 
            success = False
    if len(uniqueToRef) > 0:
            print("=== Difference: Some labels only appear in ref file:")
            from pprint import pprint 
            pprint(uniqueToRef) 
            success = False

    # Step 2: Numerical comparison
    nComp=0
    goodComp=0
    for label,results in refResult.items():
        isOK = True
        valRef = results[0]
        obsRef = results[1]
        try:
            valTest = testResult[label][0]
            obsTest = testResult[label][1]
        except KeyError:
             continue # should already have been printed in the check above 
        epsilon = 1.e-6 
        if abs(valTest-valRef) > epsilon:
             print ( f"=== Difference: Label {label} has value {valTest} in test, and {valRef} in ref - Diff is {valTest - valRef}")
             isOK=False 
        if abs(obsTest-obsRef) > 0:
             print ( f"=== Difference: Label {label} has counts {obsTest} in test, and {obsRef} in ref - Diff is {obsTest - obsRef}")
             isOK=False 
        success = success and isOK
        if isOK:
             goodComp+=1
        nComp+=1
        if lines > 0 and nComp >= lines: 
             break 
    print (f" Compared {nComp} entries, of which {goodComp} agree.")
    return success 

def compareResiduals(refFile,testFile,lines):
    refResult = parseResiduals(refFile,lines)
    testResult = parseResiduals(testFile,lines)
    success = True 

    # Step 1: Check for same cycles in file
    uniqueToTest = [label for label in testResult.keys() if label not in refResult.keys()]
    uniqueToRef = [label for label in refResult.keys() if label not in testResult.keys()]
    if len(uniqueToTest) > 0:
            print("=== Difference: Some fit cycles only appear in test file:")
            from pprint import pprint 
            pprint(uniqueToTest) 
            success = False
    if len(uniqueToRef) > 0:
            print("=== Difference: Some fit cycles only appear in ref file:")
            from pprint import pprint 
            pprint(uniqueToRef) 
            success = False

    # Now validate each cycle separately 
    for cycle, resRef in refResult.items():
         print (f" -> Validating residuals for cycle {cycle}") 
         success = compareResidualsCycle(resRef, testResult[cycle],lines) and success
    return success 

def compareLogs(refFile,testFile,lines):
    refResult = parseLogFile(refFile)
    testResult = parseLogFile(testFile)
    success = True
    # check for chi2 
    epsilon = 1.e-6 
    if abs(refResult["Chi2"]-testResult["Chi2"]) > epsilon:
         r = refResult["Chi2"]
         t = testResult["Chi2"]
         print(f"=== Difference: Fit Chi2 changed - from ref {r} to test {t}, diff is {r-t}")
         success = False 
    else:
         r = refResult["Chi2"]
         t = testResult["Chi2"]
         print (f" Chi2 in good agreement - diff is {r-t}")
    if len(refResult["convergence"]) != len(testResult["convergence"]):
           r = len(refResult["convergence"])
           t = len(testResult["convergence"])
           print (f"=== Difference: Number of iterations changed - from ref {r} to test {t}")
           success = False 
    nComp = 0
    goodComp = 0
    for i,results in enumerate(refResult["convergence"]):
        isOK = True
        fcnRef = results[0]
        deltaExpRef = results[1]
        nRejectRef = results[2]
        try:
            fcnTest = testResult["convergence"][i][0]
            deltaExpTest = testResult["convergence"][i][1]
            nRejectTest = testResult["convergence"][i][2]
        except KeyError:
             continue # should already have been printed in the check above 
        
        if abs(fcnTest-fcnRef) > 0:
             print ( f"=== Difference: Iteration {i} has functon value {fcnTest} in test, and {fcnRef} in ref - Diff is {fcnTest - fcnRef}")
             isOK=False 

        if abs(deltaExpTest-deltaExpRef) > epsilon:
             print ( f"=== Difference: Iteration {i} has expected delta {deltaExpTest} in test, and {deltaExpRef} in ref - Diff is {deltaExpTest - deltaExpRef}")
             isOK=False 

        if abs(nRejectTest-nRejectRef) > 0:
             print ( f"=== Difference: Iteration {i} has reject count {nRejectTest} in test, and {nRejectRef} in ref - Diff is {nRejectTest - nRejectRef}")
             isOK=False 


        success = success and isOK
        if isOK:
             goodComp+=1
        nComp+=1
        if lines > 0 and nComp >= lines: 
             break 
    print (f" Compared {nComp} fit iterations, of which {goodComp} agree.")
    return success 

def compareResidualsCycle(refResult,testResult,lines):
    success = True 

    # Step 1: Check for same labels in file
    uniqueToTest = [label for label in testResult.keys() if label not in refResult.keys()]
    uniqueToRef = [label for label in refResult.keys() if label not in testResult.keys()]
    if len(uniqueToTest) > 0:
            print("=== Difference: Some labels only appear in test file:")
            from pprint import pprint 
            pprint(uniqueToTest) 
            success = False
    if len(uniqueToRef) > 0:
            print("=== Difference: Some labels only appear in ref file:")
            from pprint import pprint 
            pprint(uniqueToRef) 
            success = False

    # Step 2: Numerical comparison
    epsilon = 1.e-4 
    nComp=0
    goodComp=0
    for label,results in refResult.items():
        isOK = True
        obsRef = results[0]
        medianRef = results[1]
        rmsRef = results[2]
        sigmaRef = results[3]
        try:
            obsTest = results[0]
            medianTest = results[1]
            rmsTest = results[2]
            sigmaTest = results[3]
        except KeyError:
             continue # should already have been printed in the check above 
        
        if abs(obsTest-obsRef) > 0:
             print ( f"=== Difference: Label {label} has counts {obsTest} in test, and {obsRef} in ref - Diff is {obsTest - obsRef}")
             isOK=False 

        if abs(medianTest-medianRef) > epsilon:
             print ( f"=== Difference: Label {label} has median {medianTest} in test, and {medianRef} in ref - Diff is {medianTest - medianRef}")
             isOK=False 

        if abs(rmsTest-rmsRef) > epsilon:
             print ( f"=== Difference: Label {label} has rms {rmsTest} in test, and {rmsRef} in ref - Diff is {rmsTest - rmsRef}")
             isOK=False 

        if abs(sigmaTest-sigmaRef) > epsilon:
             print ( f"=== Difference: Label {label} has sigma {sigmaTest} in test, and {sigmaRef} in ref - Diff is {sigmaTest - sigmaRef}")
             isOK=False 

        success = success and isOK
        if isOK:
             goodComp+=1
        nComp+=1
        if lines > 0 and nComp >= lines: 
             break 
    print (f" Compared {nComp} entries, of which {goodComp} agree.")
    return success 

if __name__ == "__main__":
    from os import path 
    import sys
    parser = getParser()
    args = parser.parse_args()

    if args.mode is None:
        runMode = path.splitext(args.ref)[-1][1:]
        print (f"No run mode selected - guessing '{runMode}' from file extension")
    else:
        runMode = args.mode

    success = True
    if runMode == "res":
        success = compareResults(args.ref,args.test,args.lines)
    elif runMode == "mon":
        success = compareResiduals(args.ref,args.test,args.lines)
    elif runMode == "log":
        success = compareLogs(args.ref,args.test,args.lines)
    
    if not success:
         print ("!!!!!! Ref and Test show differences. !!!!!!")
         sys.exit(99) # used in the CI to detect a succesful runs with differences
    else:
         print (" ====== Ref and Test agree. ======")
         sys.exit(0)
    import sys 
