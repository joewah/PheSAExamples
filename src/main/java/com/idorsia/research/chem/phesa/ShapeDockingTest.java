package com.idorsia.research.chem.phesa;

import com.actelion.research.chem.Molecule;
import com.actelion.research.chem.MolfileCreator;
import com.actelion.research.chem.MolfileParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.transformation.TransformationSequence;
import com.actelion.research.chem.conf.Conformer;
import com.actelion.research.chem.conf.SymmetryCorrectedRMSDCalculator;
import com.actelion.research.chem.docking.DockingFailedException;
import com.actelion.research.chem.docking.receptorpharmacophore.NegativeReceptorImageCreator;
import com.actelion.research.chem.docking.shape.ShapeDocking;
import com.actelion.research.chem.io.Mol2FileParser;
import com.actelion.research.chem.phesa.ShapeVolume;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

public class ShapeDockingTest {
	
	
	public static void dock() throws Exception {
		Mol2FileParser m2fp = new Mol2FileParser();
		MolfileParser mfp = new MolfileParser();
		StereoMolecule lig = mfp.getCompactMolecule(new File(ShapeDockingTest.class.getClassLoader().getResource("astex_1n46_ligand.mol").toURI()));
		StereoMolecule rec = m2fp.load(new File(ShapeDockingTest.class.getClassLoader().getResource("astex_1n46_protein.mol2").toURI()));
		selfDockShapeOnly(rec,lig);
		
	}
	
	public static void selfDockShapeOnly(StereoMolecule receptor, StereoMolecule nativeLigandPose) throws DockingFailedException  {
		double rmsd = -1.0;

		
		StereoMolecule toDock = new StereoMolecule(nativeLigandPose);
		toDock.ensureHelperArrays(Molecule.cHelperParities);

		TransformationSequence transform = new TransformationSequence();
		ShapeVolume bsVolume = NegativeReceptorImageCreator.create(nativeLigandPose, receptor,transform);
		ShapeDocking shapeDocking = new ShapeDocking(bsVolume,transform);

		StereoMolecule docked = shapeDocking.dock(toDock).get(0);
		
		SymmetryCorrectedRMSDCalculator rmsdCalc = new SymmetryCorrectedRMSDCalculator(new Conformer(nativeLigandPose),
					new Conformer(shapeDocking.dock(toDock).get(0)));

		rmsd = rmsdCalc.calculate();

		System.out.println(rmsd);
			
		String delim = System.getProperty("line.separator");
		File outfile = new File(receptor.getName()+ "_docked.sdf");
		   FileWriter fileWriter;
		try {
			
		    fileWriter = new FileWriter(outfile);
		    MolfileCreator mfc = new MolfileCreator(docked);
		    mfc.writeMolfile(fileWriter);
		    fileWriter.write(delim);
		    fileWriter.write("$$$$");
		    fileWriter.write(delim);
		    fileWriter.flush();
			fileWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
		e.printStackTrace();
		}
			

		

	}


}
