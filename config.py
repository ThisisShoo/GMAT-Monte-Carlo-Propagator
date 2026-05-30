script_folder = 'scripts'
data_folder = 'data'

gmat_date_format = '%d %b %Y %H:%M:%S.%f'
horizon_date_format = '%Y-%b-%d %H:%M:%S.%f'

mu_earth = 3.986004415e5  # km^3/s^2

default_params = [
    "UTCGregorian", "Earth.C3Energy",  "Earth.Altitude", "Earth.SMA",
     "Earth.ECC", "EarthMJ2000Eq.INC", "EarthMJ2000Eq.RAAN", "EarthMJ2000Eq.AOP",
     "Earth.TA", 'Sun.SMA', 'Sun.ECC', 'Heliocentric.INC', 'Heliocentric.RAAN',
     'Heliocentric.AOP', 'Sun.TA'
]

def report_config(sc_name: str, params: list = default_params,
                  filename: str = "ReportFile1.txt"):
    """Prints the configuration of the GMAT instance."""
    filename = "Report" + filename
    params = [f"{sc_name}.{param}" for param in params]
    params = ', '.join(params)
    params = "{" + params + "}"

    report_template = fr"""
%----------------------------------------
%---------- Subscribers
%----------------------------------------

Create ReportFile {filename};
GMAT {filename}.Add = {params};
GMAT {filename}.SolverIterations = Current;
GMAT {filename}.UpperLeft = [ 0.01352941176470588 0.0286493860845839 ];
GMAT {filename}.Size = [ 0.5976470588235294 0.7939972714870396 ];
GMAT {filename}.RelativeZOrder = 422;
GMAT {filename}.Maximized = false;
GMAT {filename}.Filename = '{filename}.txt';
GMAT {filename}.Precision = 16;
GMAT {filename}.WriteHeaders = true;
GMAT {filename}.LeftJustify = On;
GMAT {filename}.ZeroFill = Off;
GMAT {filename}.FixedWidth = true;
GMAT {filename}.Delimiter = ',';
GMAT {filename}.ColumnWidth = 23;
GMAT {filename}.WriteReport = true;
    """

    return report_template

if __name__ == '__main__':
    report_config("Ast2024PT5")