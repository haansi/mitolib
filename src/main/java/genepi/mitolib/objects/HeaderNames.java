package genepi.mitolib.objects;

public enum HeaderNames {
   SampleId("SampleID"),
   Position("Pos"),
   Reference("Ref"),
   VariantBase("Variant"),
   VariantLevel("Variant-Level"),
   Coverage("Coverage-Total");
   private String ColName;

   HeaderNames(String colname) {
       this.ColName = colname;
   }

   public String colname() {
       return ColName;
   }

}
