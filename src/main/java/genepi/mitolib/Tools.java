package genepi.mitolib;

import java.lang.reflect.InvocationTargetException;

import genepi.base.Toolbox;
import genepi.mitolib.bam.BAMReader;
import genepi.mitolib.contChecker.ContaminatonChecker;
import genepi.mitolib.contChecker.HaploChecker;
import genepi.mitolib.contChecker.HaploChecker2;
import genepi.mitolib.contChecker.HaploChecker2Hsd;
import genepi.mitolib.contChecker.HaploChecker3;
import genepi.mitolib.haplogrep.VariantsToHsd;
import genepi.mitolib.lofreq.LoFreqReader;
import genepi.mitolib.splitter.HeteroplasmySplitter;
import genepi.mitolib.splitter.HeteroplasmySplitterRaw;


public class Tools extends Toolbox {

	public Tools(String command, String[] args) {
		
		super(command, args);
		
	}
	
	public static void main (String[] args){
		
		if (args.length == 0) {
	        System.out.println("No arguments provided");
	      }else if (args.length > 0) {
	    	Tools tools = new Tools("java -jar mitolib.jar", args);
			
			tools.addTool("splitter", HeteroplasmySplitter.class);
			tools.addTool("splitterRaw", HeteroplasmySplitterRaw.class);
			tools.addTool("contChecker", ContaminatonChecker.class);
			tools.addTool("lofreq", LoFreqReader.class);
			tools.addTool("bam2var",BAMReader.class);
			tools.addTool("haplochecker",HaploChecker.class);
			tools.addTool("haplochecker2",HaploChecker2.class);
			tools.addTool("haplochecker2hsd",HaploChecker2Hsd.class);
			tools.addTool("haplochecker3",HaploChecker3.class);
			tools.addTool("variants2hsd",VariantsToHsd.class);
			try {
				tools.start();
			} catch (InstantiationException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (SecurityException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (NoSuchMethodException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IllegalArgumentException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (InvocationTargetException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	    }
	
	}

}
