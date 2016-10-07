package genepi.objects;


import java.util.Vector;


public class HSDEntry {
	
	public HSDEntry(String iD, String rANGE, StringBuffer pROFILES) {
		super();
		ID = iD;
		RANGE = rANGE;
		PROFILES = pROFILES;
	}
	
	public HSDEntry() {
		super();
		ID = "";
		RANGE = "";
		PROFILES = new StringBuffer();
	}

	String ID;
	String RANGE;
	StringBuffer PROFILES;
	
	public String getID() {
		return ID;
	}
	public void setID(String iD) {
		ID = iD;
	}
	public String getRANGE() {
		return RANGE;
	}
	public void setRANGE(String rANGE) {
		RANGE = rANGE;
	}
	public StringBuffer getPROFILES() {
		return PROFILES;
	}
	public void setPROFILES(StringBuffer pROFILES) {
		PROFILES = pROFILES;
	}
	
	public void appendPROFILES(String profile){
		PROFILES.append(profile+"\t");
	}
	
	public String getString(){
		return (ID+"\t"+RANGE+"\t\t"+PROFILES);
	}
}