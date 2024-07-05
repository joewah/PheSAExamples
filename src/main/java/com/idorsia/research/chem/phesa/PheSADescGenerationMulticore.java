package com.idorsia.research.chem.phesa;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.stream.Stream;

import com.actelion.research.chem.MolfileCreator;
import com.actelion.research.chem.SmilesParser;
import com.actelion.research.chem.StereoMolecule;
import com.actelion.research.chem.io.DWARFileParser;
import com.actelion.research.chem.phesa.DescriptorHandlerShape;
import com.actelion.research.chem.phesa.DescriptorHandlerShapeOneConf;
import com.actelion.research.chem.phesa.PheSAMolecule;

/**
 * Multithreaded descriptor generation: A File of SMILES is streamed into a
 * pipeline of 1) parsing SMILES into StereMolecules 2) Descriptor Calculation (involves conformer generation)
 * 3) writing the descriptors into a file
 * usage of Stream.parallel() leads to multihreaded data processing
 * to not slow down the writing, this step is done by a single thread fetching result from a queue 
 */
public class PheSADescGenerationMulticore {
	
	static AtomicInteger count = new AtomicInteger();
	private BlockingQueue<PheSAMolecule> queue;
	
	public void genDescs() throws URISyntaxException, IOException {


		Stream<String> smilesStream = Files.lines(Path.of(PheSATest.class.getClassLoader().getResource("chembl_random.txt").toURI()));

		DescriptorWriter writer = new DescriptorWriter("descriptors.txt");
		DescriptorGenerator generator = new DescriptorGenerator();
		smilesStream.parallel().map(e -> { //parallel parsing in the first and parallel descriptor generation in 2nd map step
			StereoMolecule mol = new StereoMolecule();
			try {
				new SmilesParser().parse(mol, e);
			} catch (Exception e2) {
				// TODO Auto-generated catch block
				e2.printStackTrace();
				return null;
			}
			return mol;
		}).filter(Objects::nonNull).map(e -> generator.apply(e)).forEach(e -> writer.accept(e));
		
		writer.finalizeWrite();
		

		
		
	}
	
	/**
	 * writer class that buffers incoming descriptors into a queue and writes them to disc
	 * this step is executed in a single thread
	 */
	public static class DescriptorWriter implements Consumer<PheSAMolecule> {
	        
	        private BufferedWriter writer;
	        private DescriptorHandlerShape dhs;
	        private BlockingQueue<PheSAMolecule> queue;
	        private volatile boolean running;

	        public DescriptorWriter(String filePath) {
	            dhs = new DescriptorHandlerShape();
	            queue = new LinkedBlockingQueue<>();
	            running = true;
	            try {
	                writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(filePath), "utf-8"));
	            } catch (UnsupportedEncodingException | FileNotFoundException e) {
	                e.printStackTrace();
	            }
	            new Thread(this::writeFromQueue).start();
	        }

	        private void writeFromQueue() {
	            try {
	                while (running || !queue.isEmpty()) {
	                    PheSAMolecule molecule = queue.poll();
	                    if (molecule != null) {
	                        writer.write(dhs.encode(molecule));
	                        writer.write(System.lineSeparator());
	                    }
	                }
	            } catch (IOException e) {
	                e.printStackTrace();
	            } finally {
	                try {
	                    writer.close();
	                } catch (IOException e) {
	                    e.printStackTrace();
	                }
	            }
	        }

	        public void finalizeWrite() {
	            running = false;
	        }

	        @Override
	        public void accept(PheSAMolecule t) {
	            queue.add(t);
	        }
	    }

	/**
	 * Helper class to generate descriptors from StereoMolecules, used by the parallel stream in the mapping step
	 * of the pipeline to map StereoMolecules to PheSAMolecules (PheSA descriptors)	
	 */
	public static class DescriptorGenerator implements Function<StereoMolecule, PheSAMolecule> {
		
		private DescriptorHandlerShape dhs;
	
		public DescriptorGenerator() {
			//50 conformers
			dhs = new DescriptorHandlerShape(50, 0.5);
		}

		@Override
		public PheSAMolecule apply(StereoMolecule t) {
			PheSAMolecule desc=  dhs.createDescriptor(t);
			return desc;
			
		}
		
	}

}
