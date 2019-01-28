from itertools import chain, starmap, cycle, compress

class TilesGenerator():
    
    def __init__(self, k, M, size, rank):
        self.k = k  # dimension of the contingency matrix
        self.M = M  # number of tiles across one dimension
        self.size = size  # apropo mpi
        self.rank = rank  # apropo mpi
        self._final_gen = self._checkered_tiles_gen()
        
    def __iter__(self):
        return self
    
    def __next__(self):
        try:
            return next(self._final_gen)
        except StopIteration:
            raise

    def _tiles_gen(self, sd=1):

        def embedd(*indeces):
            for i in range(indeces[-1] + sd, self.M):
                yield (*indeces, i)

        points = ((idx,) for idx in range(self.M))
        for _ in range(self.k - 1):
            points = chain.from_iterable(starmap(embedd, points))

        return points
    
    def _rank_selector(self):
        rank = self.rank
        size = self.size
        return cycle( (idx == rank for idx in range(size)) )
    
    def _checkered_tiles_gen(self):
        return compress(self._tiles_gen(), self._rank_selector())
