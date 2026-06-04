from pipeline.lines_fit import run as get_ews
from pipeline.lines_load import load_linelist
from inference.sampler import StellarInference


def run(star_id="solar"):

    # 1. get observed EWs
    ew_table = get_ews(star_id)

    ew_obs = ew_table["ew_mA"].values
    ew_err = ew_table["ew_err_mA"].values

    # 2. load line list
    line_list = load_linelist()

    # 3. run inference
    model = StellarInference(ew_obs, ew_err, line_list)

    results = model.run(nlive=500)

    return results


if __name__ == "__main__":
    run()
