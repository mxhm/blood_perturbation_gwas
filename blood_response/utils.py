from pathlib import Path
import zstandard
import lz4.frame
from pickle import dump, load, HIGHEST_PROTOCOL
import io
import pandas as pd
from tqdm.auto import tqdm
from joblib import Parallel
import time


def force_path(location, require_exists=True):
    if not isinstance(location, Path):
        location = Path(location)
    if require_exists and not location.exists():
        raise ValueError("Can't open location: {}".format(location))
    return location


def read_csv_zst(fi):
    with open(fi, 'rb') as fh:
        dctx = zstandard.ZstdDecompressor()
        with dctx.stream_reader(fh) as stream_reader:
            text_stream = io.TextIOWrapper(stream_reader, encoding='utf-8')
            table = pd.read_csv(text_stream, sep='\t')
    return table


def write_csv_zst(df, fi):
    with open(fi, 'wb') as fh:
        dctx = zstandard.ZstdCompressor()
        with dctx.stream_writer(fh) as stream_writer:
            text_stream = io.TextIOWrapper(stream_writer, encoding='utf-8')
            df.to_csv(text_stream, sep='\t', index=False)
            text_stream.flush()
            stream_writer.flush()


def lz4_load(filename):
    with lz4.frame.open(filename, mode='rb') as f:
        data = load(f)
    return data


def text_progessbar(seq, total=None):
    step = 1
    tick = time.time()
    while True:
        time_diff = time.time()-tick
        avg_speed = time_diff/step
        total_str = 'of %n' % total if total else ''
        print('step', step, '%.2f' % time_diff, 'avg: %.2f iter/sec' % avg_speed, total_str)
        step += 1
        yield next(seq)


def ParallelExecutor(use_bar='tqdm', **joblib_args):
    all_bar_funcs = {
        'tqdm': lambda args: lambda x: tqdm(x, **args, leave=False),
        'txt': lambda args: lambda x: text_progessbar(x, **args),
        'False': lambda args: iter,
        'None': lambda args: iter,
    }

    def aprun(bar=use_bar, **tq_args):
        def tmp(op_iter):
            if str(bar) in all_bar_funcs.keys():
                bar_func = all_bar_funcs[str(bar)](tq_args)
            else:
                raise ValueError("Value %s not supported as bar type" %bar)
            return Parallel(**joblib_args)(bar_func(op_iter))
        return tmp
    return aprun
