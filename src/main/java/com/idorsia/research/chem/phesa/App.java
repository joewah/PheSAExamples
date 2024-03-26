package com.idorsia.research.chem.phesa;

public class App {

	public static void main(String[] args) {
		try {
			ShapeDockingTest.dock();
			PheSATest.align();
			PheSAFlexTest.align();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
