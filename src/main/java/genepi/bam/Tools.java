package genepi.bam;

import java.lang.reflect.InvocationTargetException;

import genepi.base.Toolbox;

public class Tools extends Toolbox {

	public Tools(String command, String[] args) {
		
		super(command, args);
		
	}
	
	public static void main (String[] args){
		
		Tools tools = new Tools("java -jar greenVC.jar", args);
		
		tools.addTool("bam2var", BAMReader.class);
		
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
