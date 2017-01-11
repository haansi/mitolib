package genepi.bam;

import java.lang.reflect.InvocationTargetException;

import javax.swing.SwingUtilities;

import genepi.base.Toolbox;
import gui.greenGUI;

public class Tools extends Toolbox {

	public Tools(String command, String[] args) {
		
		super(command, args);
		
	}
	
	public static void main (String[] args){
		
		if (args.length == 0) {
	        System.out.println("No arguments - starting Gui. Start with ? for all console options");
	        SwingUtilities.invokeLater(new Runnable() {
				@Override
				public void run() {
					new greenGUI().setVisible(true);
				}
			});
	      }else if (args.length > 0) {
	    	Tools tools = new Tools("java -jar greenVC.jar", args);
			
			tools.addTool("bam2var", BAMReader.class);
			tools.addTool("haplocheck", HaploCheckReader.class);
			tools.addTool("haplocheck-mtDNA-Server", HeteroplasmyReader.class);
			tools.addTool("lofreq", LoFreqReader.class);
			
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
