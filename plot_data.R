plot_data=function(data,given_title){
  'Plots the full simulated data'
  gg=ggplot(data,aes(x=time,y=variables,color=variable))+ labs(color='Compartments') +
    geom_line(aes(y=as.numeric(E),col='Exposed')) +
    geom_line(aes(y=as.numeric(P),col='Pre-Symptomatic'))+
    geom_line(aes(y=as.numeric(I),col='Infective'))+
    geom_line(aes(y=as.numeric(A),col='Asymptomatic'))+
    geom_line(aes(y=as.numeric(J),col='Isolated'))+
    geom_line(aes(y=as.numeric(Ev),col='V-Exposed'))+
    geom_line(aes(y=as.numeric(Pv),col='V-Pre-Symptomatic'))+
    geom_line(aes(y=as.numeric(Iv),col='V-Infective'))+
    geom_line(aes(y=as.numeric(Jv),col='V-Isolated'))+
    geom_line(aes(y=as.numeric(Av),col='V-Asymptomatic'))+
    ggtitle(given_title)+
    geom_line(aes(y=as.numeric(D),col='Dead'))+
    geom_line(aes(y=as.numeric(R),col='Recovered'))+
    geom_line(aes(y=as.numeric(Tv),col='Treated'))+
    geom_line(aes(y=as.numeric(V),col='Vaccinated'))+
    geom_line(aes(y=as.numeric(S),col='Susceptible'))+
    theme(legend.key.size = unit(5, 'mm'), plot.title = element_text(size = unit(11,'cm')),
        axis.text = element_text(size = unit(8,'cm')),
        axis.title = element_text(size = unit(11,'cm')),
        legend.text = element_text(size = unit(9, 'mm')),
        legend.title =element_text(size = unit(9, 'mm')),
        legend.position = 'right')+xlab('Number of days since the 01/01/2021')+
    ylab('Size of each compartment')+geom_vline(xintercept = 29, col=1)
  print(gg)
}
