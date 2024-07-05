package com.idorsia.research.chem.phesa;

/**
 * entrypoint to run the test calculations. Tasks that should be omitted are commented out
 */
public class App {

	public static void main(String[] args) {
		try {
			//ShapeDockingTest.dock();
			//PheSATest.align();
			//PheSAFlexTest.align();
			//new PheSADescGenerationMulticore().genDescs();
			new PheSAScreenMulticore().screenLib();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
