package genepi.mitolib.objects;

public class ContaminationEntry {
String SampleId;
String majorHG;
String minorHG;
int majorFound;
int minorFound;
int majorNotFound;
int minorNotFound;
public int majorRemaining;
public int minorRemaining;
double majorFoundMean;
double minorFoundMean;
double contScore;
int contType;

public String getSampleId() {
	return SampleId;
}
public void setSampleId(String sampleId) {
	SampleId = sampleId;
}
public String getMajorId() {
	return majorHG;
}
public void setMajorId(String majorId) {
	this.majorHG = majorId;
}
public String getMinorId() {
	return minorHG;
}
public void setMinorId(String minorId) {
	this.minorHG = minorId;
}
public int getMajorFound() {
	return majorFound;
}
public void setMajorFound(int majorFound) {
	this.majorFound = majorFound;
}
public int getMinorFound() {
	return minorFound;
}
public void setMinorFound(int minorFound) {
	this.minorFound = minorFound;
}
public int getMajorNotFound() {
	return majorNotFound;
}
public void setMajorNotFound(int majorNotFound) {
	this.majorNotFound = majorNotFound;
}
public int getMinorNotFound() {
	return minorNotFound;
}
public void setMinorNotFound(int minorNotFound) {
	this.minorNotFound = minorNotFound;
}
public int getMajorRemaining() {
	return majorRemaining;
}
public void setMajorRemaining(int majorRemaining) {
	this.majorRemaining = majorRemaining;
}
public int getMinorRemaining() {
	return minorRemaining;
}
public void setMinorRemaining(int minorRemaining) {
	this.minorRemaining = minorRemaining;
}
public double getMajorFoundMean() {
	return majorFoundMean;
}
public void setMajorFoundMean(double majorFoundMean) {
	this.majorFoundMean = majorFoundMean;
}
public double getMinorFoundMean() {
	return minorFoundMean;
}
public void setMinorFoundMean(double minorFoundMean) {
	this.minorFoundMean = minorFoundMean;
}
public double getContScore() {
	return contScore;
}
public void setContScore(double contScore) {
	this.contScore = contScore;
}
public int getContType() {
	return contType;
}
public void setContType(int contType) {
	this.contType = contType;
}

	
}
