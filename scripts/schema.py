from enum import Enum, unique

from pydantic import BaseModel, Field, computed_field, field_validator


class BlastConfig(BaseModel):
    frac_identity: float = Field(
        description="Min BLAST percent identity (>0.0 and <=1.0). Default 0.90.", default=0.90
    )
    frac_aligned: float = Field(
        description="Min BLAST percent aligned (>0.0 and <=1.0). Default 0.90", default=0.90
    )
    word_size: int = Field(description="BLAST word size (5-31). Default 11", default=11)
    margin: int = Field(
        description="Margin in nucleotides for deduplicating hits by location. "
        "Hits that overlap by <= margin nucleotides are considered separate locations. Default 0.",
        default=0,
    )

    @field_validator("frac_identity", "frac_aligned")
    @classmethod
    def validate_frac(cls, value: float) -> float:
        if value <= 0.0 or value > 1.0:
            msg = f"{value} is not >=0.0 and <1.0"
            raise ValueError(msg)

        return value

    @field_validator("word_size")
    @classmethod
    def validate_word_size(cls, value: int) -> int:
        if value < 5 or value > 31:
            msg = f"{value} is not >=5 and <=31"
            raise ValueError(msg)

        return value

    @field_validator("margin")
    @classmethod
    def validate_margin(cls, value: int) -> int:
        if value < 0:
            msg = f"{value} is not >= 0"
            raise ValueError(msg)

        return value


@unique
class QueryMetaParser(Enum):
    AMRFINDER = "amr_finder"


class QueryMeta(BaseModel):
    accession: str | None = None
    locus: str | None = None
    gene_name: str | None = None
    feature_name: str | None = None
    parsed_as: QueryMetaParser | None = None


class BlastHit(BaseModel):
    # Query
    query_id: str = Field(description="name of query (from .fasta file)", alias="qseqid")
    # Subject
    subject_id: str = Field(description="name of subject (from .fasta file)", alias="sseqid")
    subject_start: int = Field(alias="sstart", description="Start of hit in subject.")
    subject_end: int = Field(alias="send", description="End of hit in subject.")
    subject_frame: int = Field(alias="sframe", description="Hit frame in subject.")
    subject_strand: str = Field(alias="sstrand", description="Hit strand in subject.")
    # Alignment
    alignment_length: int = Field(alias="length", description="Alignment length.")
    percent_identity: float = Field(
        alias="pident", description="Hit percent identity (as a fraction)."
    )
    percent_aligned: float = Field(
        alias="frac_aligned",
        description="Hit percent aligned (as a fraction). Calculated as `alignment_length / min(len(subject), len(query))`",
    )
    gaps: int = Field(alias="gaps")

    def parse_query_id_amr_finder(self) -> QueryMeta | None:
        # Early return for non-amrfinder subject id.
        if self.query_id.count("|") != 6:
            return None

        try:
            (accession, locus, _, _, gene_name, _, feature_name) = self.query_id.split("|")

            return QueryMeta(
                accession=accession,
                locus=locus,
                gene_name=gene_name,
                feature_name=feature_name,
                parsed_as=QueryMetaParser.AMRFINDER,
            )
        except ValueError:
            return None

    @computed_field
    def query_meta(self) -> QueryMeta:
        if (amr_parsed := self.parse_query_id_amr_finder()) is not None:
            return amr_parsed

        return QueryMeta()


class BlastResult(BaseModel):
    hits: list[BlastHit]
    config: BlastConfig


@unique
class BlastCols(Enum):
    # Query
    QSEQID = "qseqid"
    QLEN = "qlen"
    QSTART = "qstart"
    QEND = "qend"
    QFRAME = "qframe"
    # Subject
    SSEQID = "sseqid"
    SLEN = "slen"
    SSTART = "sstart"
    SEND = "send"
    SFRAME = "sframe"
    SSTRAND = "sstrand"
    # Alignment
    LENGTH = "length"
    PIDENT = "pident"
    GAPS = "gaps"

    @classmethod
    def to_list(cls) -> list[str]:
        return [c.value for c in cls]

    @classmethod
    def blast_str(cls) -> str:
        return " ".join(cls.to_list())
