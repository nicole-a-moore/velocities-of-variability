## vision for a cool figure
library(tidyverse)
library(plotly)

## 95% quantile colour in terrestrial env:
## min = 0.07505786
## max = 1.506056 

## min and max change:
## min = -0.03194815 per 140 years = -0.1250417/200000 days 
## max = 0.0420622 per 140 years = 0.164627/200000 days 


n.row <- 60
n.col <- 60
simulate <- function(n.row, n.col) {
  # initiate the matrix
  prob.n <- matrix(0, nrow=n.row, ncol=n.row)
  
  y.seq <- seq(from = -0.1250417, to = 0.164627 , length.out = n.row) ## change in spectral exponent
  x.seq <- seq(from = 0.07505786, to = 1.506056 , length.out = n.col) ## initial colour
  
  xx <- dnorm(x.seq, mean=n.row/2, sd=12)
  
  ## for each initial colour in spectral colour
  for (i in 1:n.row) {
    ## simulate effect of change 
    y <- append(seq(from = -0.1250417, to = 0.164627, length.out = n.row/2),  seq(to = -0.1250417, from = 0.164627, length.out = n.row/2))
    prob.n[i,] <- y*i*1000
  }
  prob.n;
}

res <- simulate(n.row, n.col)
colnames(res) <- seq(from = -0.1250417, to = 0.164627 , length.out = n.row) 
row.names(res) <- seq(from = 0.07505786, to = 1.506056 , length.out = n.col)

df <- reshape2::melt(res)
df$freq = 1:nrow(df)

fig.nc <- plot_ly(z = ~res,
                  x = seq(from = -0.1250417, to = 0.164627 , length.out = n.row),
                  y = seq(from = 0.07505786, to = 1.506056 , length.out = n.col),
                  contours = list(
                    z = list(
                      show=TRUE,
                      usecolormap=TRUE,
                      highlightcolor="#ff0000",
                      project=list(z=TRUE)
                    ),
                    y = list(
                      show=TRUE,
                      usecolormap=TRUE, 
                      highlightcolor="#ff0000",
                      project=list(y=TRUE)
                    ),
                    x = list(
                      show=TRUE,
                      usecolormap=TRUE,
                      highlightcolor="#ff0000",
                      project=list(x=TRUE)
                    )
                    
                  )
)
fig.nc <- fig.nc %>% add_surface() %>% 
  layout(
    showlegend = F,
    scene = list(
      xaxis = list(title = 'Change in spectral exponent'), 
      yaxis = list(title = 'Intitial spectral exponent'),
      zaxis = list(title = 'Mean population persistence time') 
    )
  ) %>%
  add_trace(type = 'surface',
            colorbar=list(title='Frequency of global occurrence')) 

fig.nc


plot_ly(x = seq(from = -0.1250417, to = 0.164627 , length.out = n.row),
        y = seq(from = 0.07505786, to = 1.506056 , length.out = n.col)) %>% 
  add_trace(data = df,  y=df$Var1, x=df$Var2, z=df$value, type="mesh3d") %>%
  # add_surface(
  #   z = res %>% as.matrix(),
  #   surfacecolor = res,
  #   cauto=F,
  #   cmax=max(res),
  #   cmin=min(res)
  # ) %>%
  layout(
    showlegend = F,
    scene = list(
      xaxis = list(title = 'Change in spectral exponent'), 
      yaxis = list(title = 'Intitial spectral exponent'),
      zaxis = list(title = 'Mean population persistence time') 
    )
  ) 

