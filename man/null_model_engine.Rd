\name{null_model_engine}
\alias{null_model_engine}
\title{Run null model}
\usage{
  null_model_engine(species_data, algo = "RA3",
    metric = "Pianka", N.Reps = 1000)
}
\arguments{
  \item{species_data}{a dataframe <put some guidelines in
  here>}

  \item{algo}{the algorithm to use, must be "RA1", "RA2",
  "RA3", "RA4"}

  \item{metric}{the metric used to caluclate the null
  model: choices are "Pianka", "Czekanowski", "Pianka.var",
  "Czekanowski.var", "Pianka.skew", "Czekanowski.skew";
  default is Pianka}

  \item{N.reps}{the number of replicates to run the null
  model for}
}
\description{
  Create a null model object <put a better description in
  here>
}
\examples{
\dontrun{
Put an example in here

}
}

