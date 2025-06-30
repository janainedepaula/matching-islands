# Virulence-Islands-analysis for Gipsy2
python scripts for similarities analysis in virulence islands between a specific dataset of strains. These python scripts can be used with cytoscape analysis. 

### COMPARAÇÃO DE ILHAS DE VIRULÊNCIA (OU A QUE VC QUISER) ENTRE 2 ou + LINHAGENS ###
script: match_islands.py
Os inputs do script vão ser os arquivos Virulence_Island_*.faa gerados dentro da pasta xxxx_Islands_fisher/Amino_acids. 
Neste caso, está aplicado para ilhas de virulência, mas pode ser alterado dentro do script para as demais ilhas (só mudar os nomes dos arquivos de input de virulência para o de desejo). 
Uma vez com os inputs, cria um dicionário (faa_files). 


O script faz um alinhamento de todas as sequências fastas dentro deste .faa com todas as demais dos outros arquivos .faa (ele vai fazer tanto pros arquivos do mesmo diretório, quanto para o de outro diretório - tipo:
- compara todos os fastas presentes em IHP2050_Islands_fisher/Amino_acids/Virulence_Island_1.faa vs IHP2050_Islands_fisher/Amino_acids/Virulence_Island_2.faa;
- compara todos os fastas presentes em IHP2050_Islands_fisher/Amino_acids/Virulence_Island_1.faa vs IHP2060_Islands_fisher/Amino_acids/Virulence_Island_1.faa e etc.


Ele não vai olhar somente o cabeçalho (>), vai fazer o alinhamento justamente pq se fosse olhar só como foi anotado, teria erros.
A partir disto, identifica correspondências com alta identidade de sequência (≥90%) - consegue mudar dentro do script - e relata correspondências com baixa identidade (<90%).
Indica diferença no número de genes entre as ilhas comparadas.
Aponta ilhas sem nenhuma correspondência significativa.
Gera resumos e relatórios com os dados analisados (em porcentagem de similaridade entre olhas e % de similaridade entre duas proteínas de diretórios diferentes). 


Gera dois arquivos:
- virulence_island_comparison_results: Mostra o resumo geral + comparações entre proteínas com similaridade ≥90% de diferentes inputs;
- low_identity_matches: Apresenta as demais similaridades de proteínas (<90%).


###COMPARAÇÃO PAR A PAR DE ARQUIVOS ###
Para cada par único de arquivos .faa, o script lê as sequências (genes) dos dois arquivos e para cada gene do arquivo A, encontra o gene mais similar no arquivo B.
Calcula a identidade relativa da melhor correspondência. Se a identidade for: 
- ≥ 0.9 (90%), a correspondência é armazenada como alta identidade no 'virulence_island_comparison_results'.
- < 0.9, ela entra como baixa identidade e é armazenada em 'low_identity_matches'.

### E ILHAS SEM CORRESPONDÊNCIA (ILHAS EXCLUSIVAS)? ###
O script identifica quais arquivos .faa não tiveram nenhuma correspondência com identidade ≥ 90% com qualquer outro. 
Isso é relatado no resumo como ilhas sem correspondência.

#### biblios em python usadas ###
- Bio.SeqIO: para ler os arquivos .faa e FASTA.
- Bio.Align.PairwiseAligner: para realizar o alinhamento par a par (global) entre proteínas.
	- A pontuação é definida como: match_score = 1; mismatch_score = 0, logo a "identidade" será uma fração simples: nº de matches / comprimento da maior sequência.
- os, pathlib, collections: para manipular diretórios, caminhos e estruturas auxiliares.

** O script precisa estar em uma pasta que possui ele + as duas (ou mais) pastas de output do gipsy2 tipo 'xxxx_Islands_fisher'. 
** Demora um pouco por conta do alinhamento, recomendo usar screen.

## COMO RODAR
'python3 match_islands.py' dentro do diretório que tem as pastas de interesse.

#### parte 02
script: network_similarity.py

Ele vai procurar dentro do arquivo virulence_island_comparison_results.txt gerado pelo match islands os valores de similaridades apresentados na comparação par a par. 
Depois, ele vai classificar esses valores da seguinte forma:
ilhas que possuem > 90% de similaridade;
ilhas que possuem similaridade > 50% x <90%
ilhas que possuem similaridade < 50%. 

Vai gerar um arquivo chamado: cytoscape_virulence_network_filtered.csv. 
Dentro de cytoscape_virulence_network_filtered.csv terá 4 colunas, sendo a última mostrando à qual classe aquela comparação pertence. 
Só abrir no cytoscape depois. 
  - aqui, talvez, seja bom filtrar o arquivo para mostrar somente as de interesse, pois o cytoscape_virulence_network_filtered.csv pode ser grande e atrapalhar a visualização. 

