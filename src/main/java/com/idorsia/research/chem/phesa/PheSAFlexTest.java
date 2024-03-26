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



public class PheSAFlexTest {

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
		DescriptorHandlerShape  dhs = new DescriptorHandlerShape();
		PheSAMolecule refVol = dhsOC.createDescriptor(query);
		PheSAMolecule refVolConfGen = dhs.createDescriptor(query);
		double sim = dhs.getSimilarity(refVol, refVolConfGen);
		System.out.println(sim);
		StereoMolecule[] previousAlignment = dhs.getPreviousAlignment();
		List<StereoMolecule> alignedMols = new ArrayList<>();
		alignedMols.add(previousAlignment[0]);
		alignedMols.add(previousAlignment[1]);
		int maxComparisons = 50;
		int comparisons = 0;
		dhs.setFlexible(true); //use PheSAFlex!
		while(libParser.next()) {
			if(comparisons>maxComparisons)
				break;
			StereoMolecule libMol = libParser.getMolecule();
			PheSAMolecule libVol = dhs.createDescriptor(libMol);
			sim = dhs.getSimilarity(refVol, libVol);
			System.out.println(sim);
			//getPreviousAlignment returns the query and candidate molecule in optimal alignment, it
			//they are returned in an array with length 2
			previousAlignment = dhs.getPreviousAlignment();
			alignedMols.add(previousAlignment[1]);
			comparisons++;
		}
		
		
		try(BufferedWriter structureWriter = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("aligned_flex.sdf"), "utf-8"))){
				for(StereoMolecule optimized: alignedMols) {
					MolfileCreator mfc = new MolfileCreator(optimized,false);
					mfc.writeMolfile(structureWriter);
					structureWriter.write(" \n" );
					structureWriter.write("$$$$ \n" );
				}
				
		}
		
		
	}
}