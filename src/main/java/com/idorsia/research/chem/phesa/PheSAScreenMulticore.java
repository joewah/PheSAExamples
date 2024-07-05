package com.idorsia.research.chem.phesa;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;
import java.util.function.Function;

import com.actelion.research.chem.MolfileCreator;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.alignment3d.PheSAAlignmentOptimizer.PheSASetting;
import com.actelion.research.chem.io.DWARFileParser;
import com.actelion.research.chem.io.Mol2FileParser;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.DescriptorHandlerShapeOneConf;
import com.actelion.research.chem.phesa.PheSAMolecule;

/**
 * Showcases a parallel virtual screening pipeline: encoded descriptors are streamed into a pipeline of 
 * 1) descriptor decoding into PheSAMolecules 2) result calculation (similarity + alignment) and 3) writing 
 * the best results to a SD-file
 * using Stream.parallel() leads to multithreaded PheSA calculation
 */
public class PheSAScreenMulticore {

	private BufferedWriter writer;
	static AtomicInteger atomicInt = new AtomicInteger();
	
	public PheSAScreenMulticore() throws UnsupportedEncodingException, FileNotFoundException {
		writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("results.sdf"), "utf-8"));
	}
	
	public void screenLib() {
		double simCutoff = 0.4;

		StereoMolecule queryMol = parseLigand("bace_crystal_ligand.mol2");
		PheSAMolecule query = new DescriptorHandlerShapeOneConf().createDescriptor(queryMol);
		DescriptorHandlerShape dh = new DescriptorHandlerShape();
		//align confs of query molecule to generated conformers
		System.out.println(dh.getSimilarity(query, dh.createDescriptor(queryMol)));
		
		try {
			Utils.writeMoleculeRecord(writer,dh.getPreviousAlignment()[0]);
			writer.write("$$$$ \n" );
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String filePath = "descriptors.txt";
		PheSASetting setting = new PheSASetting();
		//turn off triangles for faster screening
		setting.setUseTriangle(false);
		PheSAWriter pheSAwriter = new PheSAWriter(writer);
		PheSAScreener screener = new PheSAScreener(setting,query);
		DescriptorHandlerShape dhs = new DescriptorHandlerShape();
		try {
			Files.lines(Paths.get(filePath)).map(str -> dhs.decode(str)).parallel()
				.map(candidate -> screener.apply(candidate)).filter(e -> e.similarity>simCutoff).forEach(result -> pheSAwriter.writeResult(result));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("screened" + atomicInt.get());
		
	}
	
	public StereoMolecule parseLigand(String ligfile)  {
		StereoMolecule lig = new StereoMolecule();
		Mol2FileParser m2fp = new Mol2FileParser();
		try {
			lig = m2fp.load(new File(ShapeDockingTest.class.getClassLoader().getResource(ligfile).toURI()));
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return lig;
	}

	
	public static class PheSAWriter {
		private final String delim = System.getProperty("line.separator");
		private BufferedWriter writer;
		
		public PheSAWriter(BufferedWriter writer) {
			this.writer = writer;
		}
		
		public void writeResult(PheSAResult result)  {
			try {
				Utils.writeMoleculeRecord(writer, result.candidate);
				writeField("similarity", Double.toString(result.similarity));
				writer.write("$$$$");
				writer.write(delim);
				writer.flush();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		}
		
		   private void writeField(String name, String value) throws IOException {
			   String delim = System.getProperty("line.separator"); 
			   writer.write("> <" + name + ">");
			   writer.write(delim);
			   writer.write(value);
			   writer.write(delim);
			   writer.write(delim);
		    }
		
	}
	
	/**
	 * Screener class used in the map step to convert an input PheSADescriptor into a PheSAResult
	 * containing the similarity and the aligned molecule
	 */
	public static class PheSAScreener implements Function<PheSAMolecule,PheSAResult> {
		
		private PheSASetting setting;
		private PheSAMolecule query;
		
		public PheSAScreener(PheSASetting setting, PheSAMolecule query) {
			this.setting = setting;
			this.query = query;
		}


		@Override
		public PheSAResult apply(PheSAMolecule candidate) {
			atomicInt.getAndIncrement();
			DescriptorHandlerShape dhs = new DescriptorHandlerShape();
			dhs.setPhesaSetting(setting);
			PheSAResult result = new PheSAResult();
			try {
				result.similarity = dhs.getSimilarity(query,candidate);

				result.reference = dhs.getPreviousAlignment()[0];
				result.candidate = dhs.getPreviousAlignment()[1];
			}
			catch(Exception e) {
				System.err.println("could not align molecule");
			}
			return result;
		}
		
	}


	public static class PheSAResult {
		double similarity;
		StereoMolecule reference;
		StereoMolecule candidate;
	}
}
