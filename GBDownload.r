#####################################################
#                  GBdownload                       #
#Download DNA sequences with additional information #
#          from NCBI GenBank database               #
#                                                   #
# Camila Costa Trevisan                             #
#https://github.com/camilatrevisan/GBdownload       #
#####################################################

#GBdownload() works with two search terms, the taxon name and DNA marker name.
#As taxon you can use binomial species name, genre, family, order, etc. Please notice that accuracy with current taxonomic names depend on how the species is named on GenBank database
#GBdownload allows to retrieve up to 2000 DNA sequences in a single .fasta file, and a .csv file with information of the sampling locality (when available). For some sequences, depending on the GenBank file layout, voucher number is also retrieved.
#If the search and download sucessfully ends, you will see the following message: "Download completed with success. Check output files on your work directory."
#If the search terms return no sequences, you will see an error message, meaning that there are no sequences on GenBank matching the searck terms or there was any spelling mistake on the terms
#If you try to download more than 2000 sequences at once, you will see an error message.



#Faz a busca, baixa e junta ate 2000 sequencias em um arquivo .fasta unico e salva informacoes adicionais de localidade quando disponiveis, em um arquivo.csv
GBdownload = function(species_name, marker_name){#funcao com dois argumentos, nome da especie e nome do marcador
  ####pacotes necessarios####
  
  require(rentrez)#pacote requerido para executar entrez_search
  
  require(ape)#pacote requerido para executar read.GenBank
  
 ####Criando objeto com as especificacoes de busca####
  #transformando os argumentos da funcao gbdownload no formato aceito pela funcao  rentrez_search
  Species_Marker =  paste(#cria objeto juntando os elementos abaixo em um elemento, com a funcao paste
    species_name#nome da especie (primeiro argumento da funcao)
    ,"[Organism] AND"#parte requerida pela funcao rentrez_search
    ,marker_name#nome do marcador (segundo argumento da funcao)
    ,"[Gene]"#parte requerida pela funcao rentrez_search
    ,sep=" "#separador entre os elementos sera um espaco vazio
  )
  
  
####Executando a busca ####  
  search = entrez_search(#cria objeto com resultado da busca 
    db="nuccore"#especifica o banco de dados da NCBI desejado (GenBank)
    ,term = Species_Marker#termo de busca, objeto criado acima
    ,retmax = 2001)#limita a lista a 2001 resultados
  
  
  
  ####Criando vetores com os numeros de acesso das sequencias filtradas###
  #de 500 em 500 caso seja necessario de acordo com o resultado da busca. Foi necessario limitar a 500 por vez pois o site restringe download simultaneo de muitos dados. 
  
  if(length(search$ids) == 0) {#No caso da busca pelos termos da funcao nao encontrar resultados
    return("No sequences found under these search parameters. There are no sequences on GenBank matching the searck terms or there was any spelling mistake on the searching terms")#Encerra a funcao e retorna a mensagem de erro, de que nao ha sequencias com esses termos ou existe algum erro na grafia dos termos de busca
  }
  
  if(length(search$ids)<=500){#caso a busca retorne ate 500 resultados...
    access_numbers = na.omit(search$ids[1:500])} #...Cria um vetor com os numeros de acesso das sequencias desejadas, ate 500. Necessário o na.omit, caso contrario para qualquer que seja a quantidade de resultados da busca, access_numbers seria completado com NA ate compor 500 elementos
  
  if(length(search$ids) >= 501 & length(search$ids) <= 1000){#caso a busca retorne entre 501 e 1000 resultados...
    access_numbers = (search$ids[1:500])#Cria um vetor com os numeros de acesso das sequencias 1 a 500. Nesse caso n precisa do na.omit, pois completa as 500 sequencias
    access_numbers2 = na.omit(search$ids[501:1000])}#Cria um vetor com os numeros de acesso das sequencias desejadas de 501 a 1000, sem completar com NA o que faltar ate 1000 elementos
  
  if(length(search$ids)>=1001 & length(search$ids) <= 1500){#caso a busca retorne entre 1001 e 1500 resultados...
    access_numbers = (search$ids[1:500])#Cria um vetor com os numeros de acesso das sequencias 1 a 500. Nesse caso n precisa do na.omit, pois completa as 500 sequencias
    access_numbers2 =(search$ids[501:1000])#Cria um vetor com os numeros de acesso das sequencias 501 a 1000. Nesse caso n precisa do na.omit, pois completa as 500 sequencias
    access_numbers3 = na.omit(search$ids[1001:1500])}#Cria um vetor com os numeros de acesso das sequencias desejadas de 1001 a 1500, sem completar com NA o que faltar ate 1500
  
  if(length(search$ids)>=1501& length(search$ids) <= 2000){#caso a busca retorne entre 1501 e 2000 resultados...
    access_numbers = (search$ids[1:500])#Cria um vetor com os numeros de acesso das sequencias 1 a 500. Nesse caso n precisa do na.omit, pois completa as 500 sequencias
    access_numbers2 = (search$ids[501:1000])#Cria um vetor com os numeros de acesso das sequencias 501 a 1000. Nesse caso n precisa do na.omit, pois completa as 500 sequencias
    access_numbers3 = (search$ids[1001:1500])#Cria um vetor com os numeros de acesso das sequencias 1001 a 1500. Nesse caso n precisa do na.omit, pois completa as 500 sequencias
    access_numbers4 = na.omit(search$ids[1501:2000])}#Cria um vetor com os numeros de acesso das sequencias desejadas, sem completar com NA  
  
  if(length(search$ids)>2000){#para buscas com acima de 2000 numeros de acesso
    return("Do you really need to download all these sequences at once? Please, try to download each species or marker at a time")}#retorna mensagem sugerindo que os parametros de busca sejam revisados, para baixar menos sequencias por vez
  
  ####Baixando as sequencias de cada objeto com numeros####   
  if(length(search$ids)<=500){#para buscas que retornem ate 500 elementos
    sequences=read.GenBank(#cria objeto sequences com o resultado da funcao read.GenBank(sequencias de DNA)
      access_numbers#Executa read.GenBank no objeto que contem o numero de acesso das sequencias de 1 a 500
      ,seq.names = access.nb#mantem o numero de acesso como nome de cada sequencia
      ,species.names = FALSE#nao inclui o nome da especie ao objeto final
      ,as.character = TRUE)}#escreve as bases nitrogenadas como caracteres (a, t, g, c)
  
  
  if(length(search$ids) > 501 & length(search$ids) <= 1000){#para buscas que retornem entre 501 e 1000 elementos
    sequences=read.GenBank(#cria objeto sequences com o resultado da funcao read.GenBank(sequencias de DNA)
      access_numbers#Executa read.GenBank no objeto que contem o numero de acesso das sequencias de 1 a 500
      ,seq.names = access.nb#mantem o numero de acesso como nome de cada sequencia
      ,species.names = FALSE#nao inclui o nome da especie ao objeto final
      ,as.character = TRUE)#escreve as bases nitrogenadas como caracteres (a, t, g, c)
    sequences2=read.GenBank(#cria objeto sequences2 com o resultado da funcao read.GenBank(sequencias de DNA)
      access_numbers2#Executa read.GenBank no objeto que contem o numero de acesso das sequencias de 501 a 1000
      ,seq.names = access.nb#mantem o numero de acesso como nome de cada sequencia
      ,species.names = FALSE#nao inclui o nome da especie ao objeto final
      ,as.character = TRUE)}#escreve as bases nitrogenadas como caracteres (a, t, g, c)
  
  
  if(length(search$ids)>=1001 & length(search$ids) <= 1500){ #para buscas que retornem entre 1001 e 1500 elementos
    sequences=read.GenBank(#cria objeto sequences2 com o resultado da funcao read.GenBank(sequencias de DNA)
      access_numbers#Executa read.GenBank no objeto que contem o numero de acesso das sequencias de 1 a 500
      ,seq.names = access.nb#mantem o numero de acesso como nome de cada sequencia
      ,species.names = FALSE#nao inclui o nome da especie ao objeto final
      ,as.character = TRUE)#escreve as bases nitrogenadas como caracteres (a, t, g, c)
    
    sequences2=read.GenBank(#cria objeto sequences2 com o resultado da funcao read.GenBank(sequencias de DNA)
      access_numbers2#Executa read.GenBank no objeto que contem o numero de acesso das sequencias de 501 a 1000
      ,seq.names = access.nb#mantem o numero de acesso como nome de cada sequencia
      ,species.names = FALSE#nao inclui o nome da especie ao objeto final
      ,as.character = TRUE)#escreve as bases nitrogenadas como caracteres (a, t, g, c)
    
    sequences3=read.GenBank(#cria objeto sequences2 com o resultado da funcao read.GenBank(sequencias de DNA)
      access_numbers3#Executa read.GenBank no objeto que contem o numero de acesso das sequencias de 1001 a 1500
      ,seq.names = access.nb#mantem o numero de acesso como nome de cada sequencia
      ,species.names = FALSE#nao inclui o nome da especie ao objeto final
      ,as.character = TRUE)}#escreve as bases nitrogenadas como caracteres (a, t, g, c)
  
  
  
  if(length(search$ids)>=1501& length(search$ids) <= 2000){ #para buscas que retornem entre 1501 e 2000 elementos
    sequences=read.GenBank(#cria objeto sequences2 com o resultado da funcao read.GenBank(sequencias de DNA)
      access_numbers#Executa read.GenBank no objeto que contem o numero de acesso das sequencias de 1 a 500
      ,seq.names = access.nb#mantem o numero de acesso como nome de cada sequencia
      ,species.names = FALSE#nao inclui o nome da especie ao objeto final
      ,as.character = TRUE)#escreve as bases nitrogenadas como caracteres (a, t, g, c)
    
    sequences2=read.GenBank(#cria objeto sequences2 com o resultado da funcao read.GenBank(sequencias de DNA)
      access_numbers2#Executa read.GenBank no objeto que contem o numero de acesso das sequencias de 501 a 1000
      ,seq.names = access.nb#mantem o numero de acesso como nome de cada sequencia
      ,species.names = FALSE#nao inclui o nome da especie ao objeto final
      ,as.character = TRUE)#escreve as bases nitrogenadas como caracteres (a, t, g, c)
    
    sequences3=read.GenBank(#cria objeto sequences2 com o resultado da funcao read.GenBank(sequencias de DNA)
      access_numbers3#Executa read.GenBank no objeto que contem o numero de acesso das sequencias de 1001 a 1500
      ,seq.names = access.nb#mantem o numero de acesso como nome de cada sequencia
      ,species.names = FALSE#nao inclui o nome da especie ao objeto final
      ,as.character = TRUE)#escreve as bases nitrogenadas como caracteres (a, t, g, c)
    
    sequences4=read.GenBank(#cria objeto sequences2 com o resultado da funcao read.GenBank (sequencias de DNA)
      access_numbers4#Executa read.GenBank no objeto que contem o numero de acesso das sequencias de 1501 a 2000
      ,seq.names = access.nb#mantem o numero de acesso como nome de cada sequencia
      ,species.names = FALSE#nao inclui o nome da especie ao objeto final
      ,as.character = TRUE)}#escreve as bases nitrogenadas como caracteres (a, t, g, c)
 
  
  ####Juntando os numeros de acesso em um unico objeto####
  if(length(search$ids)<=500){#para buscas que retornem ate 500 elementos
    allsequences = sequences}#cria o objeto allsequences com as sequencias 1 a 500
  if(length(search$ids) >= 501 & length(search$ids) <= 1000){#para buscas que retornem entre 501 e 1000 elementos
    allsequences = c(sequences, sequences2)}#cria o objeto allsequences com as sequencias 1 a 500 e 501 a 100
  if(length(search$ids)>=1001 & length(search$ids) <= 1500){#para buscas que retornem entre 1001 e 1500 elementos
    allsequences = c(sequences, sequences2, sequences3)}#cria o objeto allsequences com as sequencias 1 a 500, 501 a 1000 e 100 a 1500
  if(length(search$ids)>=1501& length(search$ids) <= 2000){#para buscas que retornem entre 1501 e 2000 elementos
    allsequences = c(sequences, sequences2, sequences3, sequences4)}#cria o objeto allsequences com as sequencias 1 a 500, 501 a 1000, 1001 a 1500 e 1501 a 2000
  
  ####Salvando o arquivo com as sequencias de DNA####
  write.dna( #executa funcao que permite salvar arquivos de sequencia de DNA
    allsequences #Objeto onde estao salvas as sequencias 
    ,paste(species_name, marker_name,".fasta", sep = "_")#como padronizar o nome dos arquivos com as coisas de Species_Marker
    ,format="fasta"#formato do arquivo a ser criado com as sequencias
    ,colsep = "" #sem separacao entre os caracteres
    ,nbcol=-1)
  
 ####Acessando dados de localidade####  
 
   if(length(search$ids)<=500){#para buscas que retornem ate 500 elementos
    seq.data = entrez_summary(#cria objeto seq.data com o resultado da funcao entrez_summary, que acessa um resumo dos dados depositados nos bancos de dados do NCBI
      db="nuccore"#seleciona o banco de dados GenBank
      ,id = access_numbers)#Busca a partir dos numeros de acesso de 1 a 500
    seq.col = sapply(seq.data, function(x) {x[19]})#cria um objeto seq.col onde se aplica em todos os elementos de seq.data uma funcao que retira somente a posicao 19, que eh onde ficam os dados de localidade, quando disponiveis(ou seja, copia apenas os elementos da posicao 19 de seq.data para seq.col)
    seq.info=unlist(seq.col,  use.names = FALSE)}#seq.col eh uma lista aninhada, mais pra frente tive dificuldades em acessar o que estava nela, por isso, transformei em uma lista simples (unlist), removendo os nomes que tivessem associados aos elementos da lista (use.names = FALSE)
  
  
  
  if(length(search$ids) > 501 & length(search$ids) <= 1000){#para buscas que retornem entre 501 e 1000 elementos
    seq.data = entrez_summary(#cria objeto seq.data com o resultado da funcao entrez_summary
      db="nuccore",id = access_numbers)#seleciona o banco de dados GenBank e faz busca a partir dos numeros de acesso de 1 a 500
    seq.col = sapply(seq.data, function(x) {x[19]})#copia apenas os elementos da posicao 19 de seq.data para seq.col
    seq.info=unlist(seq.col,  use.names = FALSE)#transforma seq.col em uma lista simples seq.info
    
    seq.data2 = entrez_summary(#cria objeto seq.data2 com o resultado da funcao entrez_summary
      db="nuccore", id = access_numbers2)#seleciona o banco de dados GenBank e faz busca a partir dos numeros de acesso de 501 a 1000
    seq.col.2 = sapply(seq.data2, function(x) {x[19]})#copia apenas os elementos da posicao 19 de seq.data2 para seq.col2
    seq.info2=unlist(seq.col.2,  use.names = FALSE)}#transforma seq.col2 em uma lista simples seq.info2
  
  
  
  if(length(search$ids)>=1001 & length(search$ids) <= 1500){ #para buscas que retornem entre 1001 e 1500 elementos
    seq.data = entrez_summary(#cria objeto seq.data com o resultado da funcao entrez_summary
      db="nuccore", id = access_numbers)#seleciona o banco de dados GenBank e faz busca a partir dos numeros de acesso de 1 a 500
    seq.col = sapply(seq.data, function(x) {x[19]})#copia apenas os elementos da posicao 19 de seq.data para seq.col
    seq.info=unlist(seq.col,  use.names = FALSE)#transforma seq.col em uma lista simples seq.info
    
    seq.data2 = entrez_summary(#cria objeto seq.data2 com o resultado da funcao entrez_summary
      db="nuccore", id = access_numbers2)#seleciona o banco de dados GenBank e faz busca a partir dos numeros de acesso de 501 a 1000
    seq.col.2 = sapply(seq.data2, function(x) {x[19]})#copia apenas os elementos da posicao 19 de seq.data2 para seq.col2
    seq.info2=unlist(seq.col.2,  use.names = FALSE)#transforma seq.col2 em uma lista simples seq.info2
    
    seq.data3 = entrez_summary(#cria objeto seq.data3 com o resultado da funcao entrez_summary
      db="nuccore", id = access_numbers3)#seleciona o banco de dados GenBank e faz busca a partir dos numeros de acesso de 1001 a 1500
    seq.col.3 = sapply(seq.data3, function(x) {x[19]})#copia apenas os elementos da posicao 19 de seq.data3 para seq.col3
    seq.info3=unlist(seq.col.2,  use.names = FALSE)}#transforma seq.col3 em uma lista simples seq.info3
  
  
  
  if(length(search$ids)>=1501& length(search$ids) <= 2000){ #para buscas que retornem entre 1501 e 2000 elementos
    seq.data = entrez_summary(#cria objeto seq.data com o resultado da funcao entrez_summary
      db="nuccore", id = access_numbers)#seleciona o banco de dados GenBank e faz busca a partir dos numeros de acesso de 501 a 1000
    seq.col = sapply(seq.data, function(x) {x[19]})#copia apenas os elementos da posicao 19 de seq.data para seq.col
    seq.info=unlist(seq.col,  use.names = FALSE)#transforma seq.col em uma lista simples seq.info
    
    seq.data2 = entrez_summary(#cria objeto seq.data2 com o resultado da funcao entrez_summary
      db="nuccore", id = access_numbers2)#seleciona o banco de dados GenBank e faz busca a partir dos numeros de acesso de 501 a 1000
    seq.col.2 = sapply(seq.data2, function(x) {x[19]})#copia apenas os elementos da posicao 19 de seq.data2 para seq.col2
    seq.info2=unlist(seq.col.2,  use.names = FALSE)#transforma seq.col2 em uma lista simples seq.info2
    
    seq.data3 = entrez_summary(#cria objeto seq.data3 com o resultado da funcao entrez_summary
      db="nuccore", id = access_numbers3)#seleciona o banco de dados GenBank e faz busca a partir dos numeros de acesso de 1001 a 1500
    seq.col.3 = sapply(seq.data3, function(x) {x[19]})#copia apenas os elementos da posicao 19 de seq.data3 para seq.col3
    seq.info3=unlist(seq.col.3,  use.names = FALSE)#transforma seq.col3 em uma lista simples seq.info3
    
    seq.data4 = entrez_summary(#cria objeto seq.data4 com o resultado da funcao entrez_summary
      db="nuccore", id = access_numbers4)#seleciona o banco de dados GenBank e faz busca a partir dos numeros de acesso de 1501 a 2000
    seq.col.4 = sapply(seq.data4, function(x) {x[19]})#copia apenas os elementos da posicao 19 de seq.data4 para seq.col4
    seq.info4=unlist(seq.col.4,  use.names = FALSE)}#transforma seq.col4 em uma lista simples seq.info4
  
  
####Criando objeto com todas as informacoes de localidade juntas####  
  if(length(search$ids)<=500){#para buscas que retornem ate 500 elementos
    all.seq.info = seq.info}#cria o objeto all.seq.info com as informacoes de localidade das sequencias 1 a 500
  if(length(search$ids) >= 501 & length(search$ids) <= 1000){#para buscas que retornem entre 501 e 1000 elementos
    all.seq.info = c(seq.info, seq.info2)}#cria o objeto all.seq.info com as informacoes de localidade das sequencias 1 a 500 e 501 a 1000
  if(length(search$ids)>=1001 & length(search$ids) <= 1500){#para buscas que retornem entre 1001 e 1500 elementos
    all.seq.info = c(seq.info, seq.info2, seq.info3)}#cria o objeto all.seq.info com as informacoes de localidade das sequencias 1 a 500, 501 a 1000 e 1001 a 1500 
  if(length(search$ids)>=1501& length(search$ids) <= 2000){#para buscas que retornem entre 1501 e 2000 elementos
    all.seq.info = c(seq.info, seq.info2, seq.info3, seq.info4)}#cria o objeto all.seq.info com as informacoes de localidade das sequencias 1 a 500, 501 a 1000, 1001 a 1500 e 1501 a 2000
  
####Criando objeto com todos os numeros de acesso juntos####  
  if(length(search$ids)<=500){#para buscas que retornem ate 500 elementos
    all.access.numbers =access_numbers}#cria o objeto all.access.numbers com os numeros de acesso das sequencias 1 a 500
  if(length(search$ids) >= 501 & length(search$ids) <= 1000){#para buscas que retornem entre 501 e 1000 elementos
    all.access.numbers = c(access_numbers, access_numbers2)}#cria o objeto all.access.numbers com os numeros de acesso das sequencias 1 a 500 e 501 a 1000
  if(length(search$ids)>=1001 & length(search$ids) <= 1500){#para buscas que retornem entre 1001 e 1500 elementos
    all.access.numbers = c(access_numbers, access_numbers2, access_numbers3)}#cria o objeto all.access.numbers com os numeros de acesso das sequencias 1 a 500, 501 a 1000 e 1001 a 1500
  if(length(search$ids)>=1501& length(search$ids) <= 2000){#para buscas que retornem entre 1501 e 2000 elementos
    all.access.numbers = c(access_numbers, access_numbers2,access_numbers3, access_numbers4)} #cria o objeto all.access.numbers com os numeros de acesso das sequencias 1 a 500, 501 a 1000, 1001 a 1500 e 1501 a 2000 

  ####Salvando .csv com as informacoes de localidade####  
  SeqInfo = data.frame(#Cria um objeto SeqInfo com um dataframe 
    all.access.numbers, all.seq.info)#elementos do dataframe criado: todos os numeros de acesso e todas as informacoes de sequencias
  InfoFile = paste(species_name,"_", marker_name, ".csv", sep = "")#Organizando o  nome do arquivo a ser criado, com o nome da especie, um "_" nome do marcador, e .csv no fim, sem espacos entre esses elementos
  
  write.csv(SeqInfo, InfoFile)#Salva um arquivo .csv com o dataframe Seq.Info e as especificacoes de InfoFile como titulo
  message("Download completed with success. Check output files on your work directory.")#Ao final, exibe essa mensagem pra informar que tudo terminou bem =)
  
  
}
