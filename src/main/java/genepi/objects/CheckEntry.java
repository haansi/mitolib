package genepi.objects;


import java.util.Vector;


public class CheckEntry {
	
	 String ID;
	 int POS;
	 String REF;
  	 String BaseMajor;
	 String BaseMinor;
	 
	 double VAF;
	
	 public String getID() {
		return ID;
	}
	public void setID(String iD) {
		ID = iD;
	}
	public int getPOS() {
		return POS;
	}
	public void setPOS(int positions) {
		this.POS = positions;
	}
	public String getREF() {
		return REF;
	}
	public void setREF(String ref) {
		this.REF = ref;
	}
	public String getBaseMajor() {
		return BaseMajor;
	}
	public void setBaseMajor(String baseMajor) {
		this.BaseMajor = baseMajor;
	}
	public String getBaseMinor() {
		return BaseMinor;
	}
	public void setBaseMinor(String baseMinor) {
		this.BaseMinor = baseMinor;
	}
	public double getVAF() {
		return VAF;
	}
	public void setVAF(double vaf) {
		this.VAF = vaf;
	}

	 

	

}
