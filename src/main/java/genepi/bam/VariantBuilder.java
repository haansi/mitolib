package genepi.bam;


import genepi.util.Chromosome;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.FastaSequenceFile;

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFHeaderVersion;

import java.io.BufferedWriter;
import java.io.File;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.net.MalformedURLException;

import java.text.DecimalFormat;
import java.text.NumberFormat;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.Map.Entry;


import org.apache.commons.math3.stat.interval.AgrestiCoullInterval;
import org.apache.commons.math3.stat.interval.ClopperPearsonInterval;
import org.apache.commons.math3.stat.interval.ConfidenceInterval;
import org.apache.commons.math3.stat.interval.WilsonScoreInterval;


public class VariantBuilder {

	private String genome;
	private String reference;
	private String excludeList;
	private String outDirectory;
	private double vaf;
	private int qual;
	private boolean checkRCRS;


	public VariantBuilder(String genome) {
		this.genome = genome;
	}

	public <BaqAlt> int build() throws MalformedURLException, IOException {

		//int qualPhred = 30; // default

		String prev = "";
		checkRCRS = false;

		// create final output directory
		createFolder(outDirectory);

		final SamReader reader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT)
				.open(new File(genome));

		final SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();

		String SampleName = new File(genome).getName();

		File refFasta = new File(reference);
		FastaSequenceFile ref = new FastaSequenceFile(refFasta, true); // boolean
																		// is
																		// for
																		// truncating
																		// Names
																		// at
																		// whitespaces

		ReferenceSequence reference = ref.nextSequence();

		int size = getReferenceLength(dict, reference.length());
		String name = SampleName;
		System.out.println("NAME " + name);

		int total = 0;

		int readFWD = 0;
		int readREV = 0;
		int countAf = 0;
		int countCf = 0;
		int countGf = 0;
		int countTf = 0;
		int countDf = 0;
		int countNf = 0;
		int countAr = 0;
		int countCr = 0;
		int countGr = 0;
		int countTr = 0;
		int countDr = 0;
		int countNr = 0;
		HashMap<Integer, HashMap<String, Integer>> hm = new HashMap<>();
		Locale.setDefault(new Locale("en", "US"));
		NumberFormat formatter = new DecimalFormat("#0.0000");

		// quality threshold per base
		double het = vaf; // heteroplasmy threshold to be applied (1=100%, 0.5 =
							// 50%)
		int row = size;
		int columns = 12; // 4 base forward (A,C,G,T), N, D + all of it reverse
		int basePos =0;
		int[][] result = new int[row][columns];
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < columns; j++) {
				result[i][j] = 0;
			}
		}

		try {
			for (SAMRecord samRecord : reader) {

				// Convert read name to upper case.

				if (samRecord.getMappingQuality() > qual) {

					if (!samRecord.getReadUnmappedFlag()) {

						if (!samRecord.getDuplicateReadFlag()) {

							if (samRecord.getReadLength() > 25) {

								//Alignmentscore AS between 0 and 255 (perfect)
								
								if (!samRecord.getNotPrimaryAlignmentFlag()){
								
								if (Integer.valueOf(samRecord.getAttribute("AS") + "") > 200) {
									
									String read = samRecord.getReadString();

									for (int j = 0; j < read.length(); j++) {

										byte[] quality = samRecord.getBaseQualities();

										if (quality.length==0 || quality[j] >= qual) {
											int posBase = samRecord.getReferencePositionAtReadPosition(j + 1);
											if ((samRecord.getFlags() & 0x10) == 0x10) {

												if (checkRCRS) // if HG19 Yoruba
													posBase = checkRCRS(posBase);
												switch (read.charAt(j)) {
												case 'A':
													result[posBase % (size)][0]++;
													countAf++;
													break;
												case 'C':
													result[posBase % (size)][1]++;
													countCf++;
													break;
												case 'G':
													result[posBase % (size)][2]++;
													countGf++;
													break;
												case 'T':
													result[posBase % (size)][3]++;
													countTf++;
													break;
												case 'N':
													result[posBase % (size)][4]++;
													countNf++;
													break;
												case 'D':
													result[posBase % (size)][5]++;
													countDf++;
													break;
												default:
													break;
												}

											} else {
												if (quality.length==0 || quality[j] >= qual) {

													if (checkRCRS) // if HG19
																	// Yoruba
																	// //(reference.length()
																	// ==16571)
														posBase = checkRCRS(posBase);
													switch (read.charAt(j)) {
													case 'A':
														result[posBase % (size)][6]++;
														countAr++;
														break;
													case 'C':
														result[posBase % (size)][7]++;
														countCr++;
														break;
													case 'G':
														result[posBase % (size)][8]++;
														countGr++;
														break;
													case 'T':
														result[posBase % (size)][9]++;
														countTr++;
														break;
													case 'N':
														result[posBase % (size)][10]++;
														countNr++;
														break;
													case 'D':
														result[posBase % (size)][11]++;
														countDr++;
														break;
													default:
														break;
													}
												}
											}
										}
									}
								}

								total++;

								if ((samRecord.getFlags() & 0x10) == 16) // check
																			// forward
																			// or
																			// reverse
									readREV++;
								else
									readFWD++;

								Integer currentReferencePos = samRecord.getAlignmentStart();

								for (CigarElement cigarElement : samRecord.getCigar().getCigarElements()) {

									if (cigarElement.getOperator() == CigarOperator.D) {

										Integer cigarElementStart = currentReferencePos;
										Integer cigarElementLength = cigarElement.getLength();
										Integer cigarElementEnd = currentReferencePos + cigarElementLength;

										while (cigarElementStart < cigarElementEnd) {

											if ((samRecord.getFlags() & 0x10) == 0x10) {
												countDf++;
												result[(cigarElementEnd - cigarElementLength) % (size - 1)][11]++;
											} else {
												countDr++;
												result[(cigarElementStart - cigarElementLength + 1) % (size - 1)][5]++;
											}

											cigarElementStart += cigarElementLength;
										}

									} else
										currentReferencePos += cigarElement.getLength();

									/*
									 * TODO check insertions if
									 * (cigarElement.getOperator() ==
									 * CigarOperator.I) {
									 * 
									 * Integer cigarElementStart =
									 * currentReferencePos; Integer
									 * cigarElementLength = cigarElement
									 * .getLength();
									 * 
									 * int i = 1; while (i <=
									 * cigarElementLength) {
									 * 
									 * BasePosition basePos = counts
									 * .get(cigarElementStart + "." + i + "C");
									 * 
									 * if (basePos == null) { basePos = new
									 * BasePosition();
									 * counts.put(cigarElementStart + "." + i +
									 * "C", basePos); }
									 * 
									 * if ((samRecord.getFlags() & 0x10) ==
									 * 0x10) { basePos.addcRev(1); } else {
									 * basePos.addcFor(1); }
									 * 
									 * i++; } }
									 */
								}
							}
						}
					}
				}
			}
			}
		} catch (Exception e) {
			System.out.println("Error with BAM file");
			e.printStackTrace();
		}

		writePileup(outDirectory + File.separator + name + ".pileup", size, columns, result);
		System.out.println(outDirectory + File.separator + name);

		System.out.println("All Reads " + total + " A\tC\tG\tT\tN\tD\n" + "Forward\t" + countAf + "\t" + countCf + "\t"
				+ countGf + "\t" + countTf + "\t" + countNf + "\t" + countDf + "\n" + "Reverse\t" + countAr + "\t"
				+ countCr + "\t" + countGr + "\t" + countTr + "\t" + countNr + "\t" + countDr + "\n");

		TreeMap<Integer, String> variants = new TreeMap<>();

		for (int i = 1; i < size; i++) {
			int pos = i;
			TreeMap<Integer, String> map = new TreeMap<Integer, String>(Collections.reverseOrder());

			if (pos >= 315 && pos <= 3107)
				pos = pos - 2;
			else if (pos > 3107 && pos <= 16192)
				pos = pos - 1;
			else if (pos > 16192)
				pos = pos - 2;

			map.put(result[pos][0] + result[pos][6], "A");
			map.put(result[pos][1] + result[pos][7], "C");
			map.put(result[pos][2] + result[pos][8], "G");
			map.put(result[pos][3] + result[pos][9], "T");

			double heteroplasmy = 0.0;
			double major = 0;
			// System.out.print("\n"+i+" "+ref.charAt(i-1)+" ");

			for (Entry<Integer, String> entry : map.entrySet()) {

				if (entry.getKey() > 0) {

					if (variants.containsKey(pos)) {
						variants.put(pos, variants.get(pos) + "\t" + entry.getKey() + "\t" + entry.getValue() + "\t"
								+ formatter.format((entry.getKey() / (entry.getKey() + major))));

					} else {
						variants.put(pos, reference.getBaseString().charAt(pos - 1) + "\t" + entry.getKey() + "\t"
								+ entry.getValue());
						major = entry.getKey();
					}
				}
			}
		}

		writeVariants(outDirectory + File.separator + name + ".txt", name, formatter, het, variants);

		System.out.println(prev + " - input BAM file: " + genome);
		System.out.println(prev + " - reference size: " + reference.length());
		System.out.println(" Results saved in: " + outDirectory);

		ref.close();

		return 0;

	}

	private void writeVariants(String filename, String name, NumberFormat formatter, double het,
			TreeMap<Integer, String> variants) throws IOException {
		File fout = new File(filename);
		FileOutputStream fstream = new FileOutputStream(fout);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fstream));
		File fout1 = new File(filename + ".hsd");
		FileOutputStream fstream1 = new FileOutputStream(fout1);
		BufferedWriter bwHsd = new BufferedWriter(new OutputStreamWriter(fstream1));

		bwHsd.write(name + "\t1-16569\t\t");

		bw.write(
				"SampleID\tmtSNP\tPOS\tREFrCRS\tcountMajor\tBaseMajor\tcountMinor\tBaseMinor\tVAF\tClopperPearson_95_low\tClopperPearson_95_high\n");
		for (Entry<Integer, String> entry : variants.entrySet()) {

			String help = entry.getValue();
			StringTokenizer st = new StringTokenizer(help);
			String Ref = st.nextToken();
			Ref = Ref.toUpperCase();
			int countMajor = Integer.valueOf(st.nextToken());
			String Max = st.nextToken();

			// System.out.println(" " + entry.getValue() +" "+ entry.getKey());

			if (st.hasMoreTokens()) {
				int numberMin = Integer.valueOf(st.nextToken()); // numberMin
				String ALT = st.nextToken(); // VariantMin
				double hetlevel = Double.valueOf(st.nextToken());
				if (hetlevel >= het) {
					int n = countMajor + numberMin;
					double p = (double) (numberMin / (double) (n));

					if (Ref.equals(Max)) {
						bw.write(name + "\t" + entry.getKey() + Max + "\t" + entry.getKey() + "\t" + Ref + "\t"
								+ countMajor + "\t" + Max + "\t" + numberMin + "\t" + ALT + "\t" + hetlevel + "\t"
								+ formatter.format(getClopperPearsonInterval(n, numberMin, 0.95).getLowerBound()) + "\t"
								+ formatter.format(getClopperPearsonInterval(n, numberMin, 0.95).getUpperBound())
								+ "\n");
					} else if (!Ref.equals("N")) {
						bw.write(name + "\t" + entry.getKey() + Max + "\t" + entry.getKey() + "\t" + Ref + "\t"
								+ countMajor + "\t" + Max + "\t" + numberMin + "\t" + ALT + "\t"
								+ formatter.format(1 - hetlevel) + "\t"
								+ formatter.format(getClopperPearsonInterval(n, numberMin, 0.95).getLowerBound()) + "\t"
								+ formatter.format(getClopperPearsonInterval(n, numberMin, 0.95).getUpperBound())
								+ "\n");
						bwHsd.write(entry.getKey() + Max + "\t");
					}

				} else {
					if (!Ref.equals(Max) && hetlevel < 1 - het) {
						bw.write(name + "\t" + entry.getKey() + Max + "\t" + entry.getKey() + "\t" + Ref + "\t"
								+ countMajor + "\t" + Max + "\t" + numberMin + "\t" + ALT + "\t"
								+ formatter.format(1 - hetlevel) + "\n");
						bwHsd.write(entry.getKey() + Max + "\t");
					}
				}
			} else {
				if (!Ref.equals(Max)) {
					bw.write(name + "\t" + entry.getKey() + Max + "\t" + entry.getKey() + "\t" + Ref + "\t" + countMajor
							+ "\t" + Max + "\t" + "\t" + "\t1\n");
					bwHsd.write(entry.getKey() + Max + "\t");
				}
			}
		}
		bw.close();
		bwHsd.close();
		fstream.close();
		fstream1.close();

	}

	/**
	 * Normal Approximation Method of the Binomial Confidence Interval
	 * according: http://www.sigmazone.com/binomial_confidence_interval.htm
	 * 
	 * where p = proportion of interest n = sample size α = desired confidence
	 * z1- α/2 = “z value” for desired level of confidence z1- α/2 = 1.96 for
	 * 95% confidence z1- α/2 = 2.57 for 99% confidence z1- α/2 = 3 for 99.73%
	 * confidence
	 * 
	 * @param n
	 *            amount of observed bases
	 * @param p
	 *            ratio of heteroplasmy
	 * @param low
	 *            -1 lower confidence / 1 upper
	 * @return
	 */
	private double binConfInterval(int n, double p, int low) {
		return p + (1.96 * low) * Math.sqrt((p * (1 - p)) / n);
	}

	/**
	 * WilsonScoreInterval
	 * http://commons.apache.org/proper/commons-math/javadocs/api-3.6.1/index.html
	 * 
	 * @param coverage
	 * @param variants
	 * @param confidencelevel
	 * @return
	 */
	private ConfidenceInterval getWilsonScoreInterval(int coverage, int variants, double confidencelevel) {
		WilsonScoreInterval wi = new WilsonScoreInterval();
		ConfidenceInterval ci = wi.createInterval(coverage, variants, confidencelevel);

		return ci;
	}

	private ConfidenceInterval getAgrestiCoullInterval(int coverage, int variants, double confidencelevel) {
		AgrestiCoullInterval aci = new AgrestiCoullInterval();
		ConfidenceInterval ci = aci.createInterval(coverage, variants, confidencelevel);
		return ci;
	}

	private ConfidenceInterval getClopperPearsonInterval(int coverage, int variants, double confidencelevel) {
		ClopperPearsonInterval cpi = new ClopperPearsonInterval();
		ConfidenceInterval ci = cpi.createInterval(coverage, variants, confidencelevel);
		return ci;
	}

	private void writePileup(String filename, int size, int columns, int[][] result) throws IOException {
		File fout = new File(filename);
		FileOutputStream fstream = new FileOutputStream(fout);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(fstream));
		bw.write("POS\tA\tC\tG\tT\tN\tDEL\ta\tc\tg\tt\tn\tdel\n");

		for (int i = 0; i < size; i++) {
			bw.write((i) + "\t");
			for (int j = 0; j < columns; j++) {
				bw.write(result[i][j] + "\t");
			}
			bw.newLine();
		}
		bw.close();
		fstream.close();
	}

	private VCFHeader generateHeader(String name, String chromosome) {

		Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
		Set<String> additionalColumns = new HashSet<String>();

		headerLines.add(new VCFHeaderLine(VCFHeaderVersion.VCF4_0.getFormatString(),
				VCFHeaderVersion.VCF4_0.getVersionString()));

		headerLines.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_KEY, 1, VCFHeaderLineType.String, "Genotype"));

		additionalColumns.add(name);

		SAMSequenceDictionary sequenceDict = generateSequenceDictionary(chromosome);

		VCFHeader header = new VCFHeader(headerLines, additionalColumns);
		header.setSequenceDictionary(sequenceDict);

		return header;

	}

	private SAMSequenceDictionary generateSequenceDictionary(String chromosome) {

		SAMSequenceDictionary sequenceDict = new SAMSequenceDictionary();

		SAMSequenceRecord newSequence = new SAMSequenceRecord(chromosome, Chromosome.getChrLength(chromosome));

		sequenceDict.addSequence(newSequence);

		return sequenceDict;

	}

	private VariantContext createVC(VCFHeader header, String chrom, String rsid, List<Allele> alleles,
			List<Allele> genotype, int position) {

		final Map<String, Object> attributes = new HashMap<String, Object>();
		final GenotypesContext genotypes = GenotypesContext.create(header.getGenotypeSamples().size());

		for (final String name : header.getGenotypeSamples()) {
			final Genotype gt = new GenotypeBuilder(name, genotype).phased(false).make();
			genotypes.add(gt);
		}

		return new VariantContextBuilder("23andMe", chrom, position, position, alleles).genotypes(genotypes)
				.attributes(attributes).id(rsid).make();
	}

	public int getReferenceLength(SAMSequenceDictionary dict, int referencesize) {
		int size = 0;
		for (SAMSequenceRecord rc : dict.getSequences()) {
			if (rc.getSequenceLength() == referencesize) {
				size = referencesize + 1;
				break;
			}
		}
		return size;
	}

	private static int checkRCRS(int pos) {
		if (pos >= 315 && pos <= 3107)
			pos = pos - 2;
		else if (pos > 3107 && pos <= 16192)
			pos = pos - 1;
		else if (pos > 16192)
			pos = pos - 2;
		return pos;
	}

	public String getReference() {
		return reference;
	}

	public void setReference(String reference) {
		this.reference = reference;
	}

	public String getExcludeList() {
		return excludeList;
	}

	public void setExcludeList(String excludeList) {
		this.excludeList = excludeList;
	}

	public String getOutDirectory() {
		return outDirectory;
	}

	public void setOutDirectory(String outDirectory) {
		this.outDirectory = outDirectory;
	}

	private double getVaf() {
		return vaf;
	}

	public void setVaf(double vaf) {
		this.vaf = vaf;
	}

	private int getQual() {
		return qual;
	}

	public void setQual(int qual) {
		this.qual = qual;
	}

	/**
	 * quick copied from
	 * http://stackoverflow.com/questions/3634853/how-to-create-a-directory-in-java
	 * 
	 * @param directoryName
	 * @return
	 */
	public boolean createFolder(String directoryName) {
		File theDir = new File(directoryName);
		boolean result = false;
		// if the directory does not exist, create it
		if (!theDir.exists()) {
			System.out.println("creating directory: " + directoryName);
			result = false;

			try {
				theDir.mkdir();
				result = true;
			} catch (SecurityException se) {
				// handle it
			}
			if (result) {
				System.out.println("DIR created");
			}
		}
		return result;
	}

	
	
}
