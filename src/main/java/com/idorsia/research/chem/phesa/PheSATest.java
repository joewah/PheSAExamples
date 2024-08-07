package com.idorsia.research.chem.phesa;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.List;

import com.actelion.research.chem.MolfileCreator;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.io.DWARFileParser;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.DescriptorHandlerShapeOneConf;
import com.actelion.research.chem.phesa.PheSAMolecule;

public class PheSATest {

	/**
	 * parses a DWAR File with a 3D bioactive conformation and aligns a set of molecules on it
	 */
	public static void align() throws URISyntaxException, UnsupportedEncodingException, FileNotFoundException, IOException {
		File queryFile =  new File(PheSATest.class.getClassLoader().getResource("dude_bace1_crystal_ligand.dwar").toURI());
		File libFile = new File(PheSATest.class.getClassLoader().getResource("dude_bace1_actives_final.dwar").toURI());
		DWARFileParser queryParser = new DWARFileParser(queryFile);
		queryParser.next();
		StereoMolecule query = queryParser.getMolecule();
		DWARFileParser libParser = new DWARFileParser(libFile);
		//DescriptorHandlerOneConf should be used for the query molecule, so that the 
		//input 3D coordinates are taken and no coordinates are generated
		DescriptorHandlerShapeOneConf dhsOC = new DescriptorHandlerShapeOneConf();
		//for the candidate molecules, the standard DescriptorHandlerShape should be used, which calculates conformers
		DescriptorHandlerShape  dhs = new DescriptorHandlerShape();
		PheSAMolecule refVol = dhsOC.createDescriptor(query);
		PheSAMolecule refVolConfGen = dhs.createDescriptor(query);
		double sim = dhs.getSimilarity(refVol, refVolConfGen);
		//self-alignment check: align generated conformers of the query on the original 3D input structure
		StereoMolecule[] previousAlignment = dhs.getPreviousAlignment();
		List<StereoMolecule> alignedMols = new ArrayList<>();
		alignedMols.add(previousAlignment[0]);
		alignedMols.add(previousAlignment[1]);
		while(libParser.next()) {
			StereoMolecule libMol = libParser.getMolecule();
			PheSAMolecule libVol = dhs.createDescriptor(libMol);
			sim = dhs.getSimilarity(refVol, libVol);
			System.out.println(sim);
			//getPreviousAlignment returns the query and candidate molecule in optimal alignment, it
			//they are returned in an array with length 2
			previousAlignment = dhs.getPreviousAlignment();
			alignedMols.add(previousAlignment[1]);
		}
		
		//write the 3D alignments to an SD-File
		try(BufferedWriter structureWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("aligned.sdf"), "utf-8"))){
				for(StereoMolecule optimized: alignedMols) {
					Utils.writeMoleculeRecord(structureWriter, optimized);
					structureWriter.write("$$$$ \n" );
				}
				
		}
		
		
	}
}
