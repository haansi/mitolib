package genepi.bam;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.util.StringTokenizer;

import genepi.base.Tool;


public class LoFreqReader  extends Tool {

	public LoFreqReader(String[] args) {
		super(args);
	}
	
	
	@Override
	public void init() {

		System.out
				.println("Split mitochondrial variants according the VCF file generated with LoFreq\n\n");

	}

	@Override
	public void createParameters() {

		addParameter("in", 	"input called variants with LoFreq");
		addParameter("out", "output files for HaploGrep 2");
	}

	@Override
	public int run() {

		String in = (String) getValue("in");
		String out = (String) getValue("out");
	
		LoFreqBuilder builder = new LoFreqBuilder(in, out);

		try {
			return builder.build();
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			System.out.println("Something somewhere went terribly wrong");
			e.printStackTrace();
			return -1;
		} catch (IOException e) {
			
			System.out.println("Something somewhere went terribly wrong, again");
			// TODO Auto-generated catch block
			e.printStackTrace();
			return -1;
		}
	}

	public static void main(String[] args) {
		new BAMReader(args).start();
	}

}
