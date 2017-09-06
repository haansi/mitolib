package genepi.mitolib.lofreq;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.HashMap;

import com.google.common.io.FileWriteMode;

import genepi.base.Tool;
import genepi.mitolib.objects.CheckEntry;
import genepi.mitolib.objects.HeaderNames;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;


public class LoFreqReader extends Tool {
	


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
		double vaf =0.01;
	
		try {
			return build(in, out, vaf);
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

	public int build(String vcffile, String outfile, double vaf) throws MalformedURLException, IOException {
		
			try {
				FileWriter fw = new FileWriter(new File(outfile));
				fw.write(HeaderNames.SampleId.colname() + "\t"+ HeaderNames.Position.colname() + "\t" +HeaderNames.Reference.colname() + "\t"+
				HeaderNames.VariantBase.colname() +"\t" + HeaderNames.VariantLevel.colname() + "\t" + HeaderNames.Coverage.colname()+"\n");
				//index not required, therefore false
				VCFFileReader vcfReader = new VCFFileReader(new File(vcffile), false);
			
				CloseableIterator<VariantContext> variantIter  = vcfReader.iterator();

				try {
					while (variantIter.hasNext()) {
						final VariantContext vc = variantIter.next();
						CheckEntry entry = new CheckEntry();
						String id = vcffile;
						int pos = vc.getStart();
						String ref = vc.getReference().getBaseString();
						String alt = vc.getAlleles().get(1).getBaseString();
						double allelefreq = Double.valueOf(vc.getAttribute("AF").toString());
						String cov= (String) vc.getAttribute("DP");
						fw.write(id + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + allelefreq + "\t" + cov + "\n");
					}
				
				} catch (Exception e) {
					System.out.println("Column names not present as #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO:DP;AF;SB;DP4;CONSVAR ");
					e.printStackTrace();
				}
				fw.close();
				
		        CloserUtil.close(variantIter);
		        CloserUtil.close(vcfReader);
				
		        //genepi.mitolib.splitter.HeteroplasmySplitter.generateHSDfile(hm, outfile, vaf);
			
				System.out.println("File written");
			} catch (Exception e) {
				System.out.println("ERROR");
				e.printStackTrace();
			} 
			// Everything fine
			return 0;
		}
	
	public static void main(String args[]){
		String[] arg={"infile","outfile"};
		LoFreqReader lr = new LoFreqReader(arg);
		lr.run();
		
	}
}
