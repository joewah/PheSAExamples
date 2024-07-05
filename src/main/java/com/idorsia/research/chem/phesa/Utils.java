package com.idorsia.research.chem.phesa;

import java.io.BufferedWriter;
import java.io.IOException;

import com.actelion.research.chem.MolfileCreator;
import com.actelion.research.chem.StereoMolecule;

public class Utils {

	public static void writeMoleculeRecord(BufferedWriter writer, StereoMolecule mol) throws IOException {
		MolfileCreator mfc = new MolfileCreator(mol,false);
		mfc.writeMolfile(writer);
		writer.write(" \n" );
	}
}
