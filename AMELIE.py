import json
import requests

import urllib3


url = 'https://amelie.stanford.edu/api/vcf_api/'


def main():
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
    response = requests.post(
        url,
        verify=False,
        files={'vcfFile': open('/Users/liyaqi/PycharmProjects/PhenoAptAnalysis/realVCF/CS657.vcf', 'rb'),            # required
               # 'unaffectedDadVcf': open('dad1.vcf.gz', 'rb'),   # optional (omit if not available)
               # 'unaffectedMomVcf': open('mom1.vcf.gz', 'rb'),   # optional (omit if not available)
               # 'vcfFile2': open('kid2.vcf.gz', 'rb'),           # optional (omit if not available)
               # 'unaffectedDadVcf2': open('dad2.vcf.gz', 'rb'),  # optional (omit if not available)
               # 'unaffectedMomVcf2': open('mom2.vcf.gz', 'rb')
               }, # optional (omit if not available)
        data={'dominantAlfqCutoff': 0.1,
              'alfqCutoff': 0.5, # min. 0.1; max. 3.0 (percent)
              'filterByCount': False,
              'hmctCutoff': 1,
              'alctCutoff': 3,
              'patientName': 'Example patient',
              'patientSex': 'MALE', # or 'FEMALE', or None
              'onlyPassVariants': True,
              'filterRelativesOnlyHom': False,
              'phenotypes': ','.join(['HP:0008453', 'HP:0000912', 'HP:0002804', 'HP:0008422', 'HP:0000202', 'HP:0003422', 'HP:0004299'])})
    print(json.dumps(response.json(), indent=4))


if __name__ == "__main__":
    main()