package genepi.check;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.util.StringTokenizer;

import genepi.base.Tool;


public class HaplogrepReader  extends Tool {

	public HaplogrepReader(String[] args) {
		super(args);
	}
	
	
	@Override
	public void init() {

		System.out
				.println("Compare mitochondrial profiles from extended report in HaploGrep 2\n\n");

	}

	@Override
	public void createParameters() {

		addParameter("inHG2", 	"input HaploGrep2 extended file");
		addParameter("inVar", 	"input variant file");
		addParameter("out", "output file of contaminated Samples");
	}

	@Override
	public int run() {

		String inHG2 = (String) getValue("inHG2");
		String inVar = (String) getValue("inVar");
		String out = (String) getValue("out");
	

		HaplogrepBuilder builder = new HaplogrepBuilder(inHG2, inVar, out);


		try {
			return builder.build();
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return -1;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return -1;
		}
	}

	public static void main(String[] args) {
		new HaplogrepReader(args).start();
	}

}
