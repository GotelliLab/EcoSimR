\name{null_model_engine}
\alias{null_model_engine}
\title{Run null model}
\usage{
  null_model_engine(species_data, algo = "RA3",
    metric = "Pianka", n.reps = 1000, row.names = TRUE,
    random.seed = 0)
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

  \item{n.reps}{the number of replicates to run the null
  model.}

  \item{row.names}{Does your dataframe have row names? If
  yes, they are stripped, otherwise FALSE for data that has
  no row names}

  \item{random.seed}{Choose a seed to start your random
  number.  0 will choose a random seed, otherwise set the
  seed with any integer.}
}
\description{
  Create a null model object <put a better description in
  here>
}
\examples{
\dontrun{
## Load MacAruthur warbler data
data(macwarb)

## Run the null model
warbMod <- null_model_engine(macwarb)
## Summary and plot info
summary(warbMod)
plot(warbMod)

}
}

