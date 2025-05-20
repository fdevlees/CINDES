debug=0
from scipy import stats
import numpy as np

def print_stats(y, pred):
    from scipy.stats import pearsonr
    sq_err = (y - pred)**2

    r = pearsonr(y, pred)
    mean_err = (np.mean(sq_err), np.var(sq_err))
    perc = tuple(np.percentile(sq_err, quantile) for quantile in [25, 50, 75])

    #print "\tPearson's R: ", r
    #print "\tMean error: ", mean_err
    #print "\tPercentiles: ", perc

    ret = [ r[0], r[1], mean_err[0] ]
    ret.extend(perc)
    return ret

def order_score(order1, order2):
    r1 = np.asarray(order1)
    r2 = np.asarray(order2)
    o1 = np.argsort(r1)
    o2 = np.argsort(r2)
    print("o1:", o1)
    print("o2:", o2)
    scorefunction = lambda i,j: ( abs(i-j) )**1.0
    difs =  scorefunction( o1, o2 )
    print("difs:", difs)
    score = sum(difs)
    return score

def order_score2(order1, order2):
    r1 = np.asarray(order1)
    r2 = np.asarray(order2)
    o1 = np.argsort(r1)
    o1s= r1[o1] # sorted o1
    o2s= r2[o1] # o2 sorted in the way of o1
    o2sr=np.argsort(o2s)
    print("o2sr", o2sr)
    scorefunction = lambda i,j: ( abs(i-j) )**1.0
    rangej = np.arange(len(order1))
    difs =  scorefunction( rangej, o2sr )
    score = sum(difs)
    return score

def order_score3(order1, order2):
    return stats.spearmanr( order1, order2)[0]

def order_score4(order1, order2):
    return stats.kendalltau( order1, order2)[0]

def order_score5(order1, order2):
    ''' find the position of the best scoring structure
    here the order of the two orders is important!!!
    i assume order1 = real order and order2 = predicted order

    example:
       realscores: ( 0.3, 0.4, 0.1, 0.2 ) --> rank: [ 2 3 0 1 ]
       predscores: ( 0.5, 0.2, 0.3, 0.1 ) --> rank: [ 3 1 2 0 ]
       So: the third calculation turns out to be the best of the real values. position where real_rank is 0.
       ... now the prediction of that value has rank2. so it is two positions off from the real order.
       ... and this is accessed by: rank_pred[ rank_real.index(0) ]
    '''
    def rank(order, reverse=False):
        seq = sorted(order, reverse=reverse)
        ran = [ seq.index(v) for v in order ]
        return ran
    rank_real = rank(order1)
    rank_pred = rank(order2)
    if debug:
        print("rank_real:", rank_real)
        print("rank_pred:", rank_pred)
    offset    = rank_pred[ rank_real.index(0) ]
    return offset

if __name__ == "__main__":
    r1 = [ 3, 2, 1, 5, 4 ]
    r2 = [ 3, 4, 5, 2, 1 ]
    r3 = list(range(5))
    r4 = [ 0, 2, 4, 3, 5 ]
    r5 = [ 10, 20, 30, 40, 50, 60, 70 ]
    r6 = [ 20, 40, 10, 30, 70, 50, 60 ]
    r8 = [ 10, 40, 20, 30, 50, 70, 60 ]
    r7 = [ 20, 30, 40, 10, 70, 60, 50 ]
    score = order_score( r1, r2)
    print(score)
    score =order_score2( r1, r2)
    print(score)
    print("*"*20)
    score =order_score( r3, r4)
    print(score)
    score =order_score2( r3, r4)
    print(score)


    print("*"*20)
    #another score Kendalltau
    print("Kendall-tau:", stats.kendalltau(r1,r2)[0])
    print("Kendall-tau:", stats.kendalltau(r3,r4)[0])
    print("---")
    print("Kendall-tau:", stats.kendalltau(r5,r6)[0])
    # must be same as:
    print("Kendall-tau:", stats.kendalltau(r7,r8)[0])
    # and as:
    print("Kendall-tau:", stats.kendalltau(r8,r7)[0])

    #another score: based on score_orders1 and 2 and https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient
    # it is 1 - 6 * sum( [ i - j ]**2 ) / ( len * ( len**2 - 1 ) )
    print("*"*20)
    print("Spearmans r:", stats.spearmanr(r1,r2)[0])
    print("Spearmans r:", stats.spearmanr(r3,r4)[0])
    print("----")
    print("Spearmans r:", stats.spearmanr(r5,r6)[0])
    print("Spearmans r:", stats.spearmanr(r7,r8)[0])
    print("Spearmans r:", stats.spearmanr(r8,r7)[0])

    print("offset: ", order_score5(r5, r6 ))






