from setuptools import setup, find_packages

setup(
    name="pulsarsurveyscraper",
    version="0.1.0",
    description="Cache pulsar surveys and catalogs; allow API-like access",
    author="David Kaplan",
    author_email="kaplan@uwm.edu",
    url="",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "cache_survey=pulsarsurveyscraper.scripts.cache_survey:main",
            "search_surveys=pulsarsurveyscraper.scripts.search_surveys:main",
            "pulsarsurveyscraper_server=pulsarsurveyscraper.scripts.pulsarsurveyscraper_server:main",
        ],
    },
    python_requires=">=3.7",
    include_package_data=True,
    zip_safe=False,
)
