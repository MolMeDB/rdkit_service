import sys
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina
import re, os

class Conformers:
	def gen_conformers(self, mol, numConfs=100, maxAttempts=1000, pruneRmsThresh=0.1, useExpTorsionAnglePrefs=True, useBasicKnowledge=True, enforceChirality=True):
		ids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, maxAttempts=maxAttempts, pruneRmsThresh=pruneRmsThresh, useExpTorsionAnglePrefs=useExpTorsionAnglePrefs, useBasicKnowledge=useBasicKnowledge, enforceChirality=enforceChirality, numThreads=0)
		return list(ids)
		
	def calc_energy(self, mol, conformerId, minimizeIts):
		ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol), confId=conformerId)
		ff.Initialize()
		ff.CalcEnergy()
		results = {}
		if minimizeIts > 0:
			results["converged"] = ff.Minimize(maxIts=minimizeIts)
		results["energy_abs"] = ff.CalcEnergy()
		return results
		
	def cluster_conformers(self, mol, mode="RMSD", threshold=2.0):
		if mode == "TFD":
			dmat = TorsionFingerprints.GetTFDMatrix(mol)
		else:
			dmat = AllChem.GetConformerRMSMatrix(mol, prealigned=False)
		rms_clusters = Butina.ClusterData(dmat, mol.GetNumConformers(), threshold, isDistData=True, reordering=True)
		return rms_clusters
		
	def align_conformers(self, mol, clust_ids):
		rmslist = []
		AllChem.AlignMolConformers(mol, confIds=clust_ids, RMSlist=rmslist)
		return rmslist
			
	def run(self, mol, name, numConfs = 100, maxAttempts = 10000, pruneRmsThresh = 0.1, clusterMethod="RMSD", clusterThreshold=2.0, minimizeIterations=0):
		m = Chem.AddHs(mol)
		# generate the confomers
		conformerIds = self.gen_conformers(m, numConfs, maxAttempts, pruneRmsThresh, True, True, True)
		conformerPropsDict = {}
		for conformerId in conformerIds:
			# energy minimise (optional) and energy calculation
			conformerPropsDict[conformerId] = self.calc_energy(m, conformerId, minimizeIterations)
		# cluster the conformers
		rmsClusters = self.cluster_conformers(m, clusterMethod, clusterThreshold)

		result=[]

		clusterNumber = 0
		for cluster in rmsClusters:
			minEnergy = 9999999999999
			clusterNumber = clusterNumber+1
			rmsWithinCluster = self.align_conformers(m, cluster)
			conformer_energy = {}
			for conformerId in cluster:
				e = conformerPropsDict[conformerId]["energy_abs"]
				if e < minEnergy:
					minEnergy = e
				conformer_energy[conformerId] = float(e)
				props = conformerPropsDict[conformerId]
				props["cluster_no"] = clusterNumber
				props["cluster_centroid"] = cluster[0] + 1
				idx = cluster.index(conformerId)
				if idx > 0:
					props["rms_to_centroid"] = rmsWithinCluster[idx-1]
				else:
					props["rms_to_centroid"] = 0.0

			# Reorder conformers by energy (ASC)
			conformer_energy = dict(sorted(conformer_energy.items(), key=lambda item: item[1]))
			# Get record with minimum energy
			confId = list(conformer_energy.keys())[0]

			for sname in m.GetPropNames():
				m.ClearProp(sname)
			conformerProps = conformerPropsDict[confId]
			for key in conformerProps.keys():
				m.SetProp(key, str(conformerProps[key]))
			e = conformerProps["energy_abs"]
			if e:
				m.SetDoubleProp("energy_delta", e - minEnergy)

			result.append({
				"name": name,
				"clusterNo": conformerPropsDict[confId]["cluster_no"],
				"clusterCentroid": conformerPropsDict[confId]["cluster_centroid"],
				"rms_to_centroid": conformerPropsDict[confId]["rms_to_centroid"],
				"energy": conformer_energy[confId],
				"energyDelta": e - minEnergy,
				"sdf": Chem.MolToMolBlock(m, confId=confId)
			})

		temp = result
		all = {}

		for r in temp:
			all[r["clusterNo"]] = r["energy"]

		all = dict(sorted(all.items(), key=lambda x: x[1]))
		result = []

		for i in range(min(9, len(all))):
			NO = list(all.keys())[i]
			for t in temp:
				if t["clusterNo"] == NO:
					result.append(t)

		return result