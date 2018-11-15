#=========================================================================
# write tables of parameter sweep  outputs to files
#=========================================================================

def write_tables(plot_data, keys, zvals):

#=========================================================================
# loop over all plots to create
#=========================================================================

   with open('debug_tables.dat','w') as f:
      for key in keys:

#=========================================================================
# loop over all data sets
#=========================================================================
      
         f.write("\nData set: %s \n" % (key))
         
         iset   = 0
         for setname,dataset in plot_data.iteritems():

#=========================================================================
# extract x axis and plot options
#=========================================================================
            
            if iset == 0:
               f.write("%s    " % dataset['xlbl'])
               for x in dataset['x']:
                  f.write("%12.5f \t" % (x))
               f.write("\n")

            iset           = iset + 1
            x              = dataset['x']
            tag            = dataset['tag']

#=========================================================================
# if key exists in dataset, write it to file
#=========================================================================
   
            if key in dataset:
               f.write("%s \t" % (dataset['tag']))
               for y in dataset[key]:
                  f.write("%12.5f \t" % (y))
               f.write("\n")

   return None