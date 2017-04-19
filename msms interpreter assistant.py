#filename = 'LY-2016-08-10 29'
filename = 'LY-2014-06-12 11'
    

specific_components = { # a set of cpeific components in the solution
'NC5H4NMe2',
}

def mia(filename,dec=0):
    """MS/MS interpreter assistant"""
    def indexes(x,y, thres=0.3, min_dist=None):
        '''
        !!!!! based on PeakUtils https://bitbucket.org/lucashnegri/peakutils
        Peak detection routine.
    
        Finds the peaks in *y* by taking its first order difference. By using
        *thres* and *min_dist* parameters, it is possible to reduce the number of
        detected peaks. *y* must be signed.
    
        Parameters
        ----------
        x : list or ndarray
        y : list or ndarray (signed)
            1D amplitude data to search for peaks.
        thres : float between [0., 1.]
            Normalized threshold. Only the peaks with amplitude higher than the
            threshold will be detected.
        min_dist : int
            minimum x distance between each detected peak
    
        Returns
        -------
        ndarray
            Array containing the indexes of the peaks that were detected
        '''
        
        if isinstance(y, np.ndarray) and np.issubdtype(y.dtype, np.unsignedinteger):
            raise ValueError("y must be signed")
        if type(y) != np.ndarray: # converts to numpy array if not already
            y = np.asarray(y)
        thres = thres * (np.max(y) - np.min(y)) + np.min(y) # normalize threshold to y max
    
        # find the peaks by using the first order difference
        dy = np.diff(y) # generate a list of differences between data points
        peaks = np.where((np.hstack([dy, 0.]) < 0.)
                        & (np.hstack([0., dy]) > 0.)
                        & (y > thres))[0]
        
        if peaks.size > 1 and min_dist is not None: # if there are peaks and a minimum distance has been supplied
            highest = peaks[np.argsort(y[peaks])][::-1]
            rem = np.ones(y.size, dtype=bool)
            rem[peaks] = False
    
            for peak in highest:
                if not rem[peak]: # if the peak hasn't already been looked at
                    ind = x[peak]
                    l,r = max(0,np.searchsorted(x,ind-min_dist)),min(len(y)-1,np.searchsorted(x,ind+min_dist)) # find slice based on x values and min_dist
                    sl = slice(l,r) # create a slice object
                    #sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                    rem[sl] = True # set values in the slice to true
                    rem[peak] = False # set the peak to true
    
            peaks = np.arange(y.size)[~rem]
    
        return peaks
    
    def com_loss(dec=0,custom_losses=None):
        """takes a common loss dictionary and reduces the keys to the specified decimal place"""
        from _classes.common_losses import losses,stored_dec
        if dec > stored_dec:
            raise ValueError('The specified number of decimal places (%d) exceeds the number stored (%d)' %(dec,stored_dec))
        out = {}
        for key in losses: # round values and added to dictionary
            newkey = round(key,dec)
            if dec == 0:
                newkey = int(newkey)
            if out.has_key(newkey):
                out[newkey] += ', '
                out[newkey] += losses[key]
            else:
                out[newkey] = losses[key]
        if custom_losses is not None: # if supplied with a custom list of losses
            from _classes._Molecule import Molecule
            for item in custom_losses:
                mol = Molecule(item)
                key = round(mol.em,dec)
                if dec == 0:
                    key = int(key)
                if out.has_key(key):
                    out[key] += ', '
                    out[key] += item
                else:
                    out[key] = item
        return out
        
    def tabulate(diffs):
        """tabulates the data in the output"""
        
        string = '\t'
        for ind in inds:
            string += '%.1f\t' %x[ind]
        print string
        #string = ''
        for ind,row in enumerate(diffs):
            string = '%.1f\t' %round(x[inds[ind]],1)
            for col in diffs[ind]:
                string += '%.1f\t' %round(col,1)
            print string+'\n'
    
    def guess(diffs):
        """searches for common integer losses amoung the differences matrix and prints them"""
        loss = com_loss(0,specific_components) # grab dictionary of loss values and their probable representation
        print 'possible fragment assignments (from common losses):'
        for ind,peak in enumerate(diffs):
            for ind2,otherpeak in enumerate(diffs[ind]):
                val = int(round(otherpeak))
                if val > 0 and val in loss:
                    print `x[inds[ind]]`+' -> '+`x[inds[ind2]]`+':',val, loss[val]
    
    import numpy as np
    from _classes._mzML import mzML
    
    mzml = mzML(filename)
    x,y = mzml.sum_scans()
    
    # if not all peaks are being detected, decrease the last value handed to indexes
    inds = indexes(x,y,0.01,7)
    
    diffs = []
    for i in inds: # for each index
        difline = []
        for j in inds: # append the difference
            difline.append(x[i]-x[j])
        diffs.append(difline)
        
    tabulate(diffs) #tabulate differences in console
    guess(diffs) # guess at what the differences might mean
    
    annotations = {}
    top = max(y)
    for i in inds:
        annotations[str(x[i])] = [x[i],float(y[i])/float(top)*100.]
    from tome_v02 import plotms
    plotms([x,y],annotations=annotations,output='show')
    

if __name__ == '__main__':
    mia(filename)
    
    





