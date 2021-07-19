modulePreservation = function(
   multiData,
   multiColor,
   multiWeights = NULL,
   dataIsExpr = TRUE,
   networkType = "unsigned",
   corFnc = "cor",
   corOptions = "use = 'p'",
   referenceNetworks = 1,
   testNetworks = NULL,
   nPermutations = 100,
   includekMEallInSummary = FALSE,
   restrictSummaryForGeneralNetworks = TRUE,
   calculateQvalue = FALSE, 
   randomSeed = 12345, 
   maxGoldModuleSize = 1000, maxModuleSize = 1000, 
   quickCor = 1,
   ccTupletSize = 2,
   calculateCor.kIMall = FALSE, 
   calculateClusterCoeff = FALSE,
   useInterpolation = FALSE,
   checkData = TRUE,
   greyName = NULL,
   goldName = NULL,
   savePermutedStatistics = TRUE,
   loadPermutedStatistics = FALSE,
   permutedStatisticsFile = 
      if (useInterpolation) "permutedStats-intrModules.RData" else "permutedStats-actualModules.RData",
   plotInterpolation = TRUE,
   interpolationPlotFile = "modulePreservationInterpolationPlots.pdf",
   discardInvalidOutput = TRUE,
   parallelCalculation = FALSE,
   verbose = 1, indent = 0)
{
   # For now we use numeric labels for permutations.
   permGoldName = 0.1;
   permGreyName = 0;
  
   observed = .modulePreservationInternal(multiData, multiColor, multiWeights = multiWeights, dataIsExpr = dataIsExpr,
                     calculatePermutation = FALSE, networkType = networkType,
                     referenceNetworks = referenceNetworks, 
                     testNetworks = testNetworks, 
                     densityOnly = useInterpolation,
                     maxGoldModuleSize = maxGoldModuleSize,
                     maxModuleSize = maxModuleSize, 
                     corFnc = corFnc, corOptions = corOptions, 
                     quickCor = quickCor,
                     # calculateQuality = calculateQuality,
                     ccTupletSize = ccTupletSize,
                     calculateCor.kIMall = calculateCor.kIMall,
                     calculateClusterCoeff = calculateClusterCoeff,
                     checkData = FALSE, greyName = greyName, goldName = goldName,
                     verbose = verbose -3, indent = indent + 2);

   # Calculate preservation scores in permuted data.

    psLoaded = FALSE;
    nRegStats = 20;
    nFixStats = 3;
    permOut=list()      
    regModuleSizes = list();
    regModuleNames = list();      permutationsPresent = matrix(FALSE, nNets, nRefNets);
    interpolationUsed = matrix(FALSE, nNets, nRefNets);

    for(iref in 1:nRefNets)   # Loop over reference networks
    {
     ref = referenceNetworks[iref]
     permOut[[iref]] = list()      
     regModuleSizes[[iref]] = list();
     regModuleNames[[iref]] = list();
     nRefMods = length(unique(multiColor[[ref]]));

     for (tnet in 1:nNets)
     {
        # Retain only genes that are shared between the reference and test networks
        overlap=intersect(colnames(multiData[[ref]]$data),colnames(multiData[[tnet]]$data))
        loc1=match(overlap, colnames(multiData[[ref]]$data))
        loc2=match(overlap, colnames(multiData[[tnet]]$data))
        refName = paste("ref_", setNames[ref],sep="")
        colorRef = multiColor[[ref]][loc1]
       
        datRef=multiData[[ref]]$data[loc1, loc1]
        datTest=multiData[[tnet]]$data[loc2, loc2]

        testName=setNames[tnet]
        nRefGenes = ncol(datRef);

        colorTest=NA 
        
        name=paste(refName,"vs",testName,sep="")
        obsModSizes=list()
        nObsMods = rep(0, 2);
        tab = table(colorRef);
        nObsMods[1] = length(tab);
        obsModSizes[[1]]=tab[names(tab)!=greyName]

        permRefColors = colorRef;
        nPermMods = nObsMods[1];

        permExpr = multiData(datRef, datTest);
        names(permExpr) = setNames[c(ref, tnet)];
        permOut[[iref]][[tnet]]=list(
                  regStats = array(NA, dim = c(nPermMods+1, 
                                               nRegStats, nPermutations)),
                  fixStats = array(NA, dim = c(nObsMods[[1]], nFixStats, nPermutations)));

        # Perform actual permutations
        permColors = list();
        permColors[[1]] = permRefColors;
        permColors[[2]] = NA;
        permColorsForAcc = list();
        permColorsForAcc[[1]] = colorRef;
        permColorsForAcc[[2]] = colorTest;

         combineCalculations = function(...)
         {
           list(...);
         }
         seed = sample(1e8, 1);
         datout = foreach(perm = 1:nPermutations, .combine = combineCalculations,
                          .multicombine = TRUE, .maxcombine = nPermutations+10)%dopar% 
         {
                set.seed(seed + perm + perm^2); 
                .modulePreservationInternal(permExpr, permColors,
                                            dataIsExpr = dataIsExpr,
                                            calculatePermutation = TRUE,
                                            multiColorForAccuracy = permColorsForAcc,
                                            networkType = networkType,
                                            corFnc = corFnc, corOptions = corOptions,
                                            referenceNetworks = 1, 
                                            testNetworks = list(2),
                                            densityOnly = useInterpolation,
                                            maxGoldModuleSize = maxGoldModuleSize,
                                            maxModuleSize = maxModuleSize, quickCor = quickCor,
                                            ccTupletSize = ccTupletSize,
                                            calculateCor.kIMall = calculateCor.kIMall,
                                            calculateClusterCoeff = calculateClusterCoeff,
                                            # calculateQuality = calculateQuality,
                                            greyName = greyName, goldName = goldName,
                                            checkData = FALSE,
                                            verbose = verbose -3, indent = indent + 3)
         }
         for (perm in 1:nPermutations)
         {
           permOut[[iref]][[tnet]]$regStats[, , perm] = as.matrix(
                     cbind(datout[[perm]] [[1]]$quality[[2]][, -1],
                           datout[[perm]] [[1]]$intra[[2]],
                           datout[[perm]] [[1]]$inter[[2]]));
           permOut[[iref]][[tnet]]$fixStats[, , perm] = as.matrix(datout[[perm]] [[1]]$accuracy[[2]]);
         }
         datout = datout[[1]] # For the name stting procedures that follow...

        regStatNames = c(colnames(datout[[1]]$quality[[2]])[-1], colnames(datout[[1]]$intra[[2]]), 
                         colnames(datout[[1]]$inter[[2]]));
        regModuleNames[[iref]][[tnet]] = rownames(datout[[1]]$quality[[2]]);
        regModuleSizes[[iref]][[tnet]] = datout[[1]]$quality[[2]][, 1]
        fixStatNames = colnames(datout[[1]]$accuracy[[2]]);
        fixModuleNames = rownames(datout[[1]]$accuracy[[2]]);
        dimnames(permOut[[iref]][[tnet]]$regStats) = list(regModuleNames[[iref]][[tnet]], regStatNames,
                                                         spaste("Permutation.", c(1:nPermutations)));
        dimnames(permOut[[iref]][[tnet]]$fixStats) = list(fixModuleNames, fixStatNames,
                                                         spaste("Permutation.", c(1:nPermutations)));
        permutationsPresent[tnet, iref] = TRUE
    }

    if (savePermutedStatistics) 
     save(regModuleSizes, regStatNames, regModuleNames, fixStatNames, 
          fixModuleNames, permutationsPresent, interpolationUsed, permOut, file=permutedStatisticsFile)

   observedQuality = list();
   observedPreservation = list();
   observedReferenceSeparability = list();
   observedTestSeparability = list();
   observedOverlapCounts = list();
   observedOverlapPvalues = list();
   observedAccuracy = list();
   Z.quality = list();
   Z.preservation = list();
   Z.referenceSeparability = list();
   Z.testSeparability = list();
   Z.accuracy = list();
   interpolationStat = c(rep(TRUE, 14), FALSE, FALSE, FALSE, rep(TRUE, 6));
   for(iref in 1:nRefNets)
   {
     observedQuality[[iref]] = list();
     observedPreservation[[iref]] = list();
     observedReferenceSeparability[[iref]] = list();
     observedTestSeparability[[iref]] = list();
     observedOverlapCounts[[iref]] = list();
     observedOverlapPvalues[[iref]] = list();
     observedAccuracy[[iref]] = list();
     Z.quality[[iref]] = list();
     Z.preservation[[iref]] = list();
     Z.referenceSeparability[[iref]] = list();
     Z.testSeparability[[iref]] = list();
     Z.accuracy[[iref]] = list();
     for (tnet in 1:nNets)
     {
       nModules = nrow(observed[[iref]]$intra[[tnet]]);
       nQualiStats = ncol(observed[[iref]]$quality[[tnet]])-1;
       nIntraStats = ncol(observed[[iref]]$intra[[tnet]]);
       accuracy = observed[[iref]]$accuracy[[tnet]]
       inter = observed[[iref]]$inter[[tnet]];
       nInterStats = ncol(inter) + ncol(accuracy);
       a2i = match(rownames(accuracy), rownames(inter));
       accuracy2 = matrix(NA, nrow(inter), ncol(accuracy));
       accuracy2[a2i, ] = accuracy;
       rownames(accuracy2) = rownames(inter);
       colnames(accuracy2) = colnames(accuracy);
       modSizes = observed[[iref]]$quality[[tnet]][, 1]
       allObsStats = cbind(observed[[iref]]$quality[[tnet]][, -1], 
                           observed[[iref]]$intra[[tnet]], 
                           accuracy2, inter);
       sepCol = match("separability.qual", colnames(observed[[iref]]$quality[[tnet]]));
       sepCol2 = match("separability.pres", colnames(observed[[iref]]$intra[[tnet]]));
       quality = observed[[iref]]$quality[[tnet]][, -sepCol];

       rankColsQuality = c(5)
       rankColsDensity = 4;
     
       ranks = apply(-quality[, rankColsQuality, drop = FALSE], 2, rank, na.last = "keep");
       medRank = apply(as.matrix(ranks), 1, median, na.rm = TRUE);
       observedQuality[[iref]][[tnet]] = cbind(moduleSize = modSizes,
                                               medianRank.qual = medRank,
                                               quality[, -1]);
       preservation = cbind(observed[[iref]]$intra[[tnet]][, -sepCol2], inter );
    
       ranksDensity = apply(-preservation[, rankColsDensity, drop = FALSE], 2, rank, na.last = "keep");
       medRankDensity = apply(as.matrix(ranksDensity), 1, median, na.rm = TRUE);
       
       connSummaryInd = c(7,10) # in case of adjacency input only cor.kIM and cor.Adj which sits in the cor.cor slot
       
       ranksConnectivity = apply(-preservation[, connSummaryInd, drop = FALSE], 2, rank, na.last = "keep");
       medRankConnectivity = apply(as.matrix(ranksConnectivity), 1, median, na.rm = TRUE);
       medRank = apply(cbind(ranksDensity, ranksConnectivity), 1, median, na.rm = TRUE);
       observedPreservation[[iref]][[tnet]] = cbind(moduleSize = modSizes, 
                                                    medianRank.pres = medRank,
                                                    medianRankDensity.pres = medRankDensity,
                                                    medianRankConnectivity.pres = medRankConnectivity,
                                                    preservation);
       observedAccuracy[[iref]][[tnet]] = cbind(moduleSize = modSizes[a2i], accuracy);
       observedReferenceSeparability[[iref]][[tnet]] = cbind(moduleSize = modSizes, 
                                                     observed[[iref]]$quality[[tnet]][sepCol]);
       observedTestSeparability[[iref]][[tnet]] = cbind(moduleSize = modSizes, 
                                                observed[[iref]]$intra[[tnet]][sepCol2]);
       
       observedOverlapCounts[[iref]][[tnet]] = observed[[iref]]$overlapTables[[tnet]]$countTable;
       observedOverlapPvalues[[iref]][[tnet]] = observed[[iref]]$overlapTables[[tnet]]$pTable;

       nAllStats = ncol(allObsStats);
       zAll = matrix(NA, nModules, nAllStats);
       rownames(zAll) = rownames(observed[[iref]]$intra[[tnet]])
       colnames(zAll) = paste("Z.", colnames(allObsStats), sep="");
       logModSizes = log(modSizes);
       goldRowPerm = match(goldName, regModuleNames[[iref]][[tnet]]);
       goldRowObs = match(goldName, rownames(inter));
       fixInd = 1;
       regInd = 1;
       for (stat in 1:nAllStats) 
         if (interpolationStat[stat])
         {
           means = c(apply(permOut[[iref]][[tnet]]$regStats[, regInd, , drop = FALSE], 
                     c(1:2), mean, na.rm = TRUE));
           SDs = c(apply(permOut[[iref]][[tnet]]$regStats[, regInd, , drop = FALSE], 
                     c(1:2), sd, na.rm = TRUE));
           z = ( allObsStats[, stat] - means) / SDs;
           if (any(is.finite(z)))
           {
              finite = is.finite(z)
              z[finite][SDs[finite]==0] = max(abs(z[finite]), na.rm = TRUE) *
                            sign(allObsStats[, stat] - means)[SDs[finite]==0]
           }
           zAll[, stat] = z;
           regInd = regInd + 1;
         } else {
            fixInd = fixInd + 1
         }
       zAll = as.data.frame(zAll);
       sepCol = match("Z.separability.qual", colnames(zAll)[1:nQualiStats]);
       zQual = zAll[, c(1:nQualiStats)][ , -sepCol];
       summaryColsQuality = rankColsQuality -1 # quality also contains module sizes, Z does not
       summZ = apply(zQual[, summaryColsQuality, drop = FALSE], 1, median, na.rm = TRUE); 
       Z.quality[[iref]][[tnet]] = data.frame(cbind(moduleSize = modSizes, Zsummary.qual = summZ, 
                                                    zQual));
       Z.referenceSeparability[[iref]][[tnet]] = 
             data.frame(cbind(moduleSize = modSizes, zAll[sepCol]));
       st = nQualiStats + 1; en = nQualiStats + nIntraStats + nInterStats;
       sepCol2 = match("Z.separability.pres", colnames(zAll)[st:en]);
       accuracyCols = match(c("Z.accuracy", "Z.minusLogFisherP", "Z.coClustering"), colnames(zAll)[st:en]);
       zPres = zAll[, st:en][, -c(sepCol2, accuracyCols)];
       summaryColsPreservation = list(summaryColsQuality, connSummaryInd);
       nGroups = length(summaryColsPreservation)
       summZMat = matrix(0, nrow(zPres), nGroups);
       for (g in 1:nGroups)
         summZMat[, g] = apply(zPres[, summaryColsPreservation[[g]], drop = FALSE], 1, median, na.rm = TRUE);
       colnames(summZMat) = c("Zdensity.pres", "Zconnectivity.pres");
       summZ = apply(summZMat, 1, mean, na.rm = TRUE);
       Z.preservation[[iref]][[tnet]] = data.frame(cbind(moduleSize = modSizes, Zsummary.pres = summZ, 
                                                         summZMat, zPres));
       Z.testSeparability[[iref]][[tnet]] = 
             data.frame(cbind(moduleSize = modSizes, zAll[nQualiStats + sepCol2]));
       Z.accuracy[[iref]][[tnet]] = 
             data.frame(cbind(moduleSize = modSizes, zAll[st:en][accuracyCols]));
     }
     names(observedQuality[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedPreservation[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedReferenceSeparability[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedTestSeparability[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedAccuracy[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedOverlapCounts[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(observedOverlapPvalues[[iref]]) = paste("inColumnsAlsoPresentIn", sep=".",
                                            setNames);
     names(Z.quality[[iref]]) = paste("inColumnsAlsoPresentIn", sep = ".",
                                      setNames);
     names(Z.preservation[[iref]]) = paste("inColumnsAlsoPresentIn", sep = ".",
                                      setNames);
     names(Z.referenceSeparability[[iref]]) = paste("inColumnsAlsoPresentIn", sep = ".",
                                                    setNames);
     names(Z.testSeparability[[iref]]) = paste("inColumnsAlsoPresentIn", sep = ".",
                                                    setNames);
     names(Z.accuracy[[iref]]) = paste("inColumnsAlsoPresentIn", sep = ".",
                                                    setNames);
   } 

   names(observedQuality) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedPreservation) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedReferenceSeparability) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedTestSeparability) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedAccuracy) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedOverlapCounts) = paste("ref", setNames[referenceNetworks], sep=".");
   names(observedOverlapPvalues) = paste("ref", setNames[referenceNetworks], sep=".");
   names(Z.quality) = paste("ref", setNames[referenceNetworks], sep=".");
   names(Z.preservation) = paste("ref", setNames[referenceNetworks], sep=".");
   names(Z.referenceSeparability) = paste("ref", setNames[referenceNetworks], sep=".");
   names(Z.testSeparability) = paste("ref", setNames[referenceNetworks], sep=".");
   names(Z.accuracy) = paste("ref", setNames[referenceNetworks], sep=".");

   summaryIndQuality = list(c(3:6));
   summaryIndPreservation = list(c(5:8), connSummaryInd + 4, c(3,4)); #4 = 1 module size + 3 summary indices

   output=list(quality = list(observed = observedQuality, Z = Z.quality, 
                              log.p = p.quality, 
                              log.pBonf = pBonf.quality,
                              q = q.quality ),
               preservation = list(observed = observedPreservation, Z = Z.preservation,
                                   log.p = p.preservation,
                                   log.pBonf = pBonf.preservation,
                                   q = q.preservation),
               accuracy = list(observed = observedAccuracy, Z = Z.accuracy,
                               log.p = p.accuracy,
                               log.pBonf = pBonf.accuracy,
                               q = q.accuracy,
                               observedCounts = observedOverlapCounts,
                               observedFisherPvalues = observedOverlapPvalues),
               referenceSeparability = list(observed = observedReferenceSeparability,
                                            Z = Z.referenceSeparability,
                                            log.p = p.referenceSeparability,
                                            log.pBonf = pBonf.referenceSeparability,
                                            q = q.referenceSeparability),
               testSeparability = list(observed = observedTestSeparability,
                                       Z = Z.testSeparability,
                                       log.p = p.testSeparability,
                                       log.pBonf = pBonf.testSeparability,
                                       q = q.testSeparability),
               permutationDetails = list(permutedStatistics = permOut, 
                                         interpolationModuleSizes = regModuleSizes,
                                         interpolationStatNames = regStatNames,
                                         permutationsPresent = permutationsPresent,
                                         interpolationUsed = interpolationUsed)
              );
   return(output)
}


#=====================================================================================================
#
# .modulePreservationInternal
#
#=====================================================================================================

# Calculate module preservation scores for a given multi-expression data set.

.accuracyStatistics = function(colorRef, colorTest, ccTupletSize, greyName, pEpsilon)
{
   colorRefLevels = sort(unique(colorRef));
   nRefMods = length(colorRefLevels);
   refModSizes = table(colorRef);

   accuracy=matrix(NA,nRefMods ,3)
   colnames(accuracy)=c("accuracy", "minusLogFisherP", "coClustering");
   rownames(accuracy)=colorRefLevels;

   nRefGenes = length(colorRef); # also equals nTestGenes
   overlap = list(countTable = NA, pTable = NA);
   list(accuracy = accuracy, overlapTable = overlap);
}

#=================================================================================================
#
# .modulePreservationInternal
#
#=================================================================================================

# multiData contains either expression data or adjacencies; which one is indicated in dataIsExpr

.modulePreservationInternal = function(multiData, multiColor, 
                              multiWeights,
                              dataIsExpr,
                              calculatePermutation,
                              multiColorForAccuracy = NULL,
                              networkType = "signed", 
                              corFnc = "cor", 
                              corOptions = "use = 'p'",
                              referenceNetworks=1,
                              testNetworks,
                              densityOnly = FALSE,
                              maxGoldModuleSize = 1000, 
                              maxModuleSize = 1000, quickCor = 1,
                              ccTupletSize,
                              calculateCor.kIMall,
                              calculateClusterCoeff,
                              # calculateQuality = FALSE, 
                              checkData = TRUE,
                              greyName,
                              goldName,
                              pEpsilon = 1e-200,
                              verbose = 1, indent = 0)
{

   spaces = indentSpaces(indent);

   nNets = length(multiData);
   nGenes = sapply(multiData, sapply, ncol);

   nType = charmatch(networkType, .networkTypes);

   setNames = names(multiData);

   keepGenes = list();
   for (s in 1:nNets) 
      keepGenes[[s]] = rep(TRUE, nGenes[s])
        
   datout=list()	
  for (iref in 1:length(referenceNetworks))
  {
    ref = referenceNetworks[iref]
    accuracy=list()
    quality = list();
    interPres=list()
    intraPres=list()
    overlapTables = list();

    netPresent = rep(FALSE, nNets)
    for (tnet in testNetworks[[iref]])
    {
      overlap=intersect(colnames(multiData[[ref]]$data),colnames(multiData[[tnet]]$data))
      loc1=match(overlap, colnames(multiData[[ref]]$data))
      loc2=match(overlap, colnames(multiData[[tnet]]$data))
      datTest=multiData[[tnet]]$data[loc2, loc2]
      datRef=multiData[[ref]]$data[loc1, loc1]
      colorRef=multiColor[[ref]][loc1]
  
       # Accuracy measures

       if (calculatePermutation)
       {
          colorRefAcc = sample(colorRefAcc);
          colorTestAcc = sample(colorRefAcc);
       }
       x = .accuracyStatistics(colorRefAcc, colorTestAcc, 
                              ccTupletSize = ccTupletSize, greyName = greyName, pEpsilon = pEpsilon);
       accuracy[[tnet]] = x$accuracy;
       overlapTables[[tnet]] = x$overlapTable;

       # From now on we work with colorRef; colorTest is not needed anymore.

       # Restrict each module to at most maxModuleSize genes..

       colorRefLevels = sort(unique(colorRef));
       nRefMods = length(colorRefLevels);
       nRefGenes = length(colorRef);

       # Check that the gold module is not too big. In particular, the gold module must not contain all valid
       # genes, since in such a case the random sampling makes no sense for density-based statistics.

       goldModSize = maxGoldModuleSize
       if (goldModSize > nRefGenes/2) goldModSize = nRefGenes/2;

       # ..step 1: gold module. Note that because of the above the gold module size is always smaller than
       # nRefGenes.

       goldModR = sample(nRefGenes, goldModSize)
       goldModT = sample(nRefGenes, goldModSize)

       goldRef = datRef[goldModR, goldModR];
       goldRefP = datRef[goldModT, goldModT];
       goldTest = datTest[goldModT, goldModT];

       # ..step 2: proper modules and grey

       keepGenes = rep(TRUE, nRefGenes);
       for (m in 1:nRefMods)
       {
         inModule = colorRef == colorRefLevels[m]
         nInMod = sum(inModule)
         if(nInMod > maxModuleSize)
         {
            sam = sample(nInMod, maxModuleSize)
            keepGenes[inModule] = FALSE;
            keepGenes[inModule][sam] = TRUE;
         }
       }

       # Create the permuted data sets
       if (sum(keepGenes) < nRefGenes)
       {
         colorRef = colorRef[keepGenes]
         datRef = datRef[keepGenes, keepGenes]
         nRefGenes = length(colorRef);
         keepPerm = sample(nRefGenes, sum(keepGenes));
         datTest = datTest[keepPerm, keepPerm];
         datRefP = datRef[keepPerm, keepPerm];
       } else {
         perm = sample(c(1:nRefGenes));
         datRefP = datRef[perm, perm];
         datTest = datTest[perm, perm];
       }
       datRef = .combineAdj(datRef, goldRef);
       datRefP = .combineAdj(datRefP, goldRefP);
       datTest = .combineAdj(datTest, goldTest);
       if (!is.null(rownames(datRef))) rownames(datRef) = make.unique(rownames(datRef));
       if (!is.null(rownames(datRefP))) rownames(datRefP) = make.unique(rownames(datRefP));
       if (!is.null(rownames(datTest))) rownames(datTest) = make.unique(rownames(datTest));
       }
       if (!is.null(colnames(datRef))) colnames(datRef) = make.unique(colnames(datRef));
       if (!is.null(colnames(datRefP))) colnames(datRefP) = make.unique(colnames(datRefP));
       if (!is.null(colnames(datTest))) colnames(datTest) = make.unique(colnames(datTest));

       gold = rep(goldName, goldModSize)
       colorRef_2 = c(as.character(colorRef),gold)
       colorLevels = sort(unique(colorRef_2));
       opt = list(corFnc = corFnc, corOptions = corOptions, quickCor = quickCor, 
                  nType = nType, 
                  MEgold = MEgold, MEgrey = MEgrey, 
                  densityOnly = densityOnly, calculatePermutation = calculatePermutation,
                  calculateCor.kIMall = calculateCor.kIMall, 
                  calculateClusterCoeff = calculateClusterCoeff);

      stats = .coreCalcForAdj(datRef, datRefP, datTest, colorRef_2, opt);
      interPresNames = spaste(corFnc, c(".kIM", ".kME", ".kIMall", ".adj", ".clusterCoeff", ".MAR"));
      measureNames = c("propVarExplained", "meanKIM", "separability", 
                       "meanSignAwareCorDat", "meanAdj", "meanClusterCoeff", "meanMAR");

       name1=paste(setNames[[ref]],"_vs_",setNames[[tnet]],sep="")  
       quality[[tnet]] = cbind(stats$modSizes, 
                               stats$proVar[, 1], 
                               if (dataIsExpr) stats$meanSignAwareKME[, 1] else stats$meankIM[, 1],
                               stats$Separability[, 1], stats$MeanSignAwareCorDat[,1],
                               stats$MeanAdj[, 1], stats$meanClusterCoeff[, 1], 
                               stats$meanMAR[, 1]);
       intraPres[[tnet]]=cbind(stats$proVar[, 2], 
                               if (dataIsExpr) stats$meanSignAwareKME[, 2] else stats$meankIM[, 2],
                               stats$Separability[, 2], stats$MeanSignAwareCorDat[, 2], 
                               stats$MeanAdj[, 2], stats$meanClusterCoeff[, 2],
                               stats$meanMAR[, 2])
       colnames(quality[[tnet]]) = c("moduleSize", paste(measureNames, "qual", sep="."));
       rownames(quality[[tnet]]) = colorLevels
       colnames(intraPres[[tnet]]) = paste(measureNames, "pres", sep=".");
       rownames(intraPres[[tnet]]) = colorLevels
       names(intraPres)[tnet]=paste(name1,sep="")
       quality[[tnet]] = as.data.frame(quality[[tnet]]);
       intraPres[[tnet]] = as.data.frame(intraPres[[tnet]]);
       interPres[[tnet]]= as.data.frame(cbind(stats$corkIM, stats$corkME, stats$corkMEall, stats$ICORdat,
                                              stats$corCC, stats$corMAR))
       colnames(interPres[[tnet]])=interPresNames;
       rownames(interPres[[tnet]])=colorLevels
       names(interPres)[[tnet]]=paste(name1,sep="")
       netPresent[tnet] = TRUE;
   } # of for (test in testNetworks[[iref]])

   datout[[iref]]=list(netPresent = netPresent, quality = quality, 
                       intra = intraPres, inter = interPres, accuracy = accuracy,
                       overlapTables = overlapTables)
 } # of for (iref in 1:length(referenceNetworkss))
 names(datout)=setNames[referenceNetworks]
  return(datout)
        
}

.combineAdj = function(block1, block2)
{
  n1 = ncol(block1);
  n2 = ncol(block2);
  comb = matrix(0, n1+n2, n1+n2);
  comb[1:n1, 1:n1] = block1;
  comb[(n1+1):(n1+n2), (n1+1):(n1+n2)] = block2;
  try( {colnames(comb) = c(colnames(block1), colnames(block2)) }, silent = TRUE);
  comb;
}


# This function is basically copied from the file networkConcepts.R

.computeLinksInNeighbors = function(x, imatrix){x %*% imatrix %*% x}
.computeSqDiagSum = function(x, vec) { sum(x^2 * vec) };

.clusterCoeff = function(adjmat1)
{
  # diag(adjmat1)=0
  no.nodes=dim(adjmat1)[[1]]
  nolinksNeighbors <- c(rep(-666,no.nodes))
  total.edge <- c(rep(-666,no.nodes))
  maxh1=max(as.dist(adjmat1) ); minh1=min(as.dist(adjmat1) );
  nolinksNeighbors <- apply(adjmat1, 1, .computeLinksInNeighbors, imatrix=adjmat1)
  subTerm = apply(adjmat1, 1, .computeSqDiagSum, vec = diag(adjmat1));
  plainsum  <- colSums(adjmat1)
  squaresum <- colSums(adjmat1^2)
  total.edge = plainsum^2 - squaresum
  #CChelp=rep(-666, no.nodes)
  CChelp=ifelse(total.edge==0,0, (nolinksNeighbors-subTerm)/total.edge)
  CChelp
} 

# This function assumes that the diagonal of the adjacency matrix is 1
.MAR = function(adjacency)
{
  denom = apply(adjacency, 2, sum)-1;
  mar = (apply(adjacency^2, 2, sum) - 1)/denom;
  mar[denom==0] = NA;
  mar;
}

#===================================================================================================
#
# Core calculation for adjacency
#
#===================================================================================================

# A few supporting functions first:

.getSVDs = function(data, colors)
{
  colorLevels = levels(factor(colors))
  nMods =length(colorLevels)
  svds = list();
  for (m in 1:nMods)
  {
    modGenes = (colors==colorLevels[m])
    modAdj = data[modGenes, modGenes];
    if (sum(is.na(modAdj))>0)
    {
      seed = .Random.seed;
      modAdj = impute.knn(modAdj)$data;
      .Random.seed <<- seed;
    }
    svds[[m]] = svd(modAdj, nu=1, nv=0);
    svds[[m]]$u = c(svds[[m]]$u);
    if (sum(svds[[m]]$u, na.rm = TRUE) < 0) svds[[m]]$u = -svds[[m]]$u;
  }
  svds;
}

.kIM = function(adj, colors, calculateAll = TRUE)
{
  colorLevels = levels(factor(colors))
  nMods =length(colorLevels)
  nGenes = length(colors);
  kIM = matrix(NA, nGenes, nMods);
  if (calculateAll)
  {
     for (m in 1:nMods)
     {
       modGenes = colors==colorLevels[m];
       kIM[, m] = apply(adj[, modGenes, drop = FALSE], 1, sum, na.rm = TRUE);
       kIM[modGenes, m] = kIM[modGenes, m] - 1;
     }
  } else {
     for (m in 1:nMods)
     {
       modGenes = colors==colorLevels[m];
       kIM[modGenes, m] = apply(adj[modGenes, modGenes, drop = FALSE], 1, sum, na.rm = TRUE) - 1;
     }
  }
  kIM;
}


# Here is the main function

# Summary: 
#     PVE: from svd$d
#     kME: from svd$u (to make it different from kIM)
#     kIM: as usual
#     kMEall: from kIMall
#     meanSignAwareKME: from svd$u 
#     Separability: as in the paper
#     MeanSignAwareCorDat: ??
#     meanAdj: mean adjacency
#     cor.cor: replace by cor.adj

.coreCalcForAdj = function(datRef, datRefP, datTest, colors, opt)
{
  colorLevels = levels(factor(colors))
  nMods =length(colorLevels)
  nGenes = length(colors);

  svds=list()
  kIM = list();
         
  svds[[1]] = .getSVDs(datRef, colors);
  kIM[[1]] = .kIM(datRef, colors, calculateAll = opt$calculateCor.kIMall);
  svds[[2]] = .getSVDs(datRefP, colors);
  kIM[[2]] = .kIM(datRefP, colors, calculateAll = opt$calculateCor.kIMall);
  svds[[3]] = .getSVDs(datTest, colors);
  kIM[[3]] = .kIM(datTest, colors, calculateAll = opt$calculateCor.kIMall);

  proVar=matrix(NA, nMods ,2)

  modGenes = list();
  for (m in 1:nMods)
  {
    modGenes[[m]] = c(1:nGenes)[colors==colorLevels[m]];
    proVar[m, 1] = svds[[2]][[m]]$d[1]/sum(svds[[2]][[m]]$d); 
    proVar[m, 2] = svds[[3]][[m]]$d[1]/sum(svds[[3]][[m]]$d); 
  }

  corkME = rep(NA, nMods);
  corkMEall = rep(NA, nMods);
  corkIM = rep(NA, nMods);
  ICORdat = rep(NA,nMods)
   
  for(m in 1:nMods )
  {
     nModGenes=modSizes[m];
     corExpr = parse(text=paste(opt$corFnc, "(svds[[1]][[m]]$u,svds[[3]][[m]]$u",
                                             prepComma(opt$corOptions), ")"));
     corkME[m] = abs(eval(corExpr));
     corExpr = parse(text=paste(opt$corFnc, "(kIM[[1]][modGenes[[m]],m],kIM[[3]][modGenes[[m]],m] ", 
                                prepComma(opt$corOptions), ")"));
     corkIM[m] = eval(corExpr);

     adj1 = datRef[modGenes[[m]], modGenes[[m]]];
     adj2 = datTest[modGenes[[m]], modGenes[[m]]];
     corExpr = parse(text=paste(opt$corFnc, "(c(as.dist(adj1)), c(as.dist(adj2))",
                                prepComma(opt$corOptions), ")"));
     ICORdat[m] = eval(corExpr);
  }
      
  meanSignAwareKME=matrix(NA, nMods ,2)
  meankIM=matrix(NA, nMods ,2)
  for(m in 1:nMods )
  {       
     meankIM[m, 1] = mean(kIM[[2]][modGenes[[m]], m], na.rm = TRUE)
     meankIM[m, 2] = mean(kIM[[3]][modGenes[[m]], m], na.rm = TRUE)
     meanSignAwareKME[m,1]=mean(abs(svds[[2]][[m]]$u),na.rm = TRUE)
     meanSignAwareKME[m,2]=abs(mean(sign(svds[[1]][[m]]$u) * svds[[3]][[m]]$u,na.rm = TRUE))
  }
  
  MeanAdj = matrix(NA, nMods, 2);
  sepMat = array(NA, dim = c(nMods, nMods, 2));
  corCC = rep(NA, nMods);
  corMAR = rep(NA, nMods);
  meanCC = matrix(NA,nMods ,2)
  meanMAR = matrix(NA,nMods ,2)
  for (m in 1:nMods)
  {
    modAdj = datRefP[modGenes[[m]], modGenes[[m]]];
    marRefP = .MAR(modAdj);
    meanMAR[m, 1] = mean(marRefP);
    MeanAdj[m,1]=mean(as.dist(modAdj), na.rm = TRUE);

    modAdj = datRef[modGenes[[m]], modGenes[[m]]];
    marRef = .MAR(modAdj);
  
    modAdj = datTest[modGenes[[m]], modGenes[[m]]];
    marTest = .MAR(modAdj);
    MeanAdj[m,2] = mean(as.dist(modAdj), na.rm = TRUE);
     
    meanMAR[m, 2] = mean(marTest);
    corExpr = parse(text=paste(opt$corFnc, "(marRef, marTest ", prepComma(opt$corOptions), ")"));
    corMAR[m] = eval(corExpr);

    if ((m > 1) && (colorLevels[m]!=gold))
    {
       for (m2 in 1:(m-1)) if (colorLevels[m2]!=gold)
       {
          interAdj = datRefP[modGenes[[m]], modGenes[[m2]]];
          tmp = mean(interAdj, na.rm = TRUE);
          if (tmp!=0) {
            sepMat[m, m2, 1] = mean(interAdj, na.rm = TRUE)/sqrt(MeanAdj[m, 1] * MeanAdj[m2, 1]);
          } else 
            sepMat[m, m2, 1] = 0;
          sepMat[m2, m, 1] = sepMat[m, m2, 1];
          interAdj = datTest[modGenes[[m]], modGenes[[m2]]];
          tmp = mean(interAdj, na.rm = TRUE);
          if (tmp!=0) {
            sepMat[m, m2, 2] = mean(interAdj, na.rm = TRUE)/sqrt(MeanAdj[m, 2] * MeanAdj[m2, 2]);
          } else
            sepMat[m, m2, 2] = 0;
          sepMat[m2, m, 2] = sepMat[m, m2, 2];
       }
    }
  }
  Separability=matrix(NA, nMods ,2)
  notGold = colorLevels!=gold;
  for(k in 1:2)
     Separability[notGold, k]=1-apply(sepMat[notGold, notGold, k, drop = FALSE], 1, max, na.rm = TRUE)                        
  MeanSignAwareCorDat=matrix(NA,nMods ,2)
  list(modSizes = modSizes, 
       corkIM = corkIM, corkME = corkME, corkMEall = corkMEall, 
       proVar = proVar, 
       meanSignAwareKME = meanSignAwareKME,
       meankIM = meankIM,
       Separability = Separability, MeanSignAwareCorDat = MeanSignAwareCorDat, ICORdat = ICORdat,
       MeanAdj = MeanAdj, meanClusterCoeff = meanCC, meanMAR = meanMAR, corCC = corCC, corMAR = corMAR)
}
