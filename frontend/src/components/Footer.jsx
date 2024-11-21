import { Row, Col, Button, Space, Layout } from 'antd'
import githubLogo from "../assets/img/github.svg"
import folderLogo from "../assets/img/folder.svg"
import locationLogo from "../assets/img/location.svg"
import zenodoLogo from '../assets/img/zenodo.svg'
import emailLogo from "../assets/img/email.svg"

const {Footer} = Layout


export const PairpotFooter = () => {
  return (
    <Footer>
    <Row align={'stretch'}>
      <Col span={14} style={{ marginBottom: '1rem' }} offset={1}>
        <Row>
          <Col style={{ marginTop: '0.3rem' }}>
            <h3>About Pairpot</h3>
          </Col>
          <Col style={{ display: 'flex', alignItems: 'center' }}>
          <Space direction='horizontal' size={0}>
            <a href="https://github.com/lyotvincent/Pairpot">
              <img src={githubLogo} alt="GitHub Logo" style={{ height: '30px'}} />
            </a>
            <a href="src.bioxai.cn">
              <img src={folderLogo} alt="Folder Logo" style={{ height: '20px',marginRight:'5px'}} />
            </a>
            {/* <a href="http://bio.nankai.edu.cn/">
              <img src={locationLogo} alt="Location Logo" style={{ height: '1.6rem' }} />
            </a> */}
            <a href="https://zenodo.org/records/13919102">
              <img src={zenodoLogo} alt="Zenodo Logo" style={{ height: '20px', marginLeft:'-13px'}} />
            </a>
            </Space>
          </Col>
          <p style={{ fontSize: '16px' }}>
            Pairpot's data resources are available at <a href='http://src.bioxai.cn'>src</a> in bulk.
            Any issues can be reported at <a href='https://github.com/lyotvincent/Pairpot/issues'>GitHub Issue</a>.
            Any feedback will be appreciated. <b>Please contact:</b>
          </p>
          <p style={{ marginTop: '0rem', fontSize: '16px' }}>
            Zhihan Ruan{'  '} <img src={emailLogo} style={{ width: '16px', height: '16px', marginTop: 0 }} />{' '}<i>rrrzhan@mail.nankai.edu.cn</i><br />
            Jian Liu{'  '} <img src={emailLogo} style={{ width: '16px', height: '16px', marginTop: 0 }} />{' '}<i>jianliu@nankai.edu.cn</i><br />
            <img src={locationLogo} alt="Location Logo" style={{ width: '16px', height: '16px', marginTop: -5 }} />{' '}<i>State Key Laboratory of Medicinal Chemical Biology, Centre for Bioinformatics and Intelligent Medicine, College of Computer Science, Nankai University, Tianjin, China.</i>
          </p>
        </Row>
      </Col>

      <Col span={4} style={{ marginBottom: '1rem' }} offset={2}>
        <Col style={{ marginTop: '0.3rem' }}>
          <h3>Other Works</h3>
        </Col>
        <Row>
          <Button href='http://neotcrdb.bioxai.cn/home'
            style={{ backgroundColor: 'transparent', boxShadow: 'none', borderColor: 'transparent', textDecoration: 'none', marginTop: '0.5rem' }}>
            {"> NeoTCR"}
          </Button>
        </Row>
        <Row>
          <Button href='https://dupscan.sysumeg.com'
            style={{ backgroundColor: 'transparent', boxShadow: 'none', borderColor: 'transparent', textDecoration: 'none' }}>
            {"> DupScan"}
          </Button>
        </Row>
        <Row>
          <Button href='https://hic.bioaimed.com/'
            style={{ backgroundColor: 'transparent', boxShadow: 'none', borderColor: 'transparent', textDecoration: 'none' }}>
            {"> HiBrowser"}
          </Button>
        </Row>
        <Row>
          <Button href='https://www.neopeptide.cn/'
            style={{ backgroundColor: 'transparent', boxShadow: 'none', borderColor: 'transparent', textDecoration: 'none' }}>
            {"> NeoPeptide"}
          </Button>
        </Row>
      </Col>


      {/* <Col span={7} style={{ marginBottom: '1rem' }} offset={1}>
        <Col style={{ marginTop: '0.3rem' }}>
          <h3>Visitors</h3>
        </Col>
        <a href="https://www.revolvermaps.com/livestats/5v61plwc2pd/"><img src="//rf.revolvermaps.com/h/m/a/5/00ff6c/128/25/5v61plwc2pd.png" alt="Map" style={{ border: "0", width: "256px", height: "128px" }} /></a>
      </Col> */}
      <Col style={{ marginBottom: '1rem' }} span={21} offset={1}>
      <b>Citation:</b> Ruan Z., Lin F., Zhang Z., Cao J., Xiang W., Wei X. and Liu J. (2024) Pairpot: a database with real-time lasso-based analysis tailored for paired single-cell and spatial transcriptomics. Nucleic Acids Research, gkae986, &nbsp;
            <a href="https://doi.org/10.1093/nar/gkae986" target='_blank'
              style={{
                color: 'rgba(0, 0, 0, 0.6)',
                textDecoration: 'none',
                // fontSize: 24
              }}>
              https://doi.org/10.1093/nar/gkae986
            </a>
            </Col>
            <Col style={{ marginBottom: '1rem' }} span={5} offset={20}>
      <p style={{ marginTop: '-1rem', fontSize: '14px', marginLeft: '-0.3rem' }}>
        <a href="https://beian.miit.gov.cn" target='_blank'
          style={{
            color: 'rgba(0, 0, 0, 0.6)',
            textDecoration: 'none',
            fontSize: 13,
            textAlign: 'center',
          }}>
          津ICP备2024022781号-1
        </a>
      </p>
      </Col>
    </Row>
  </Footer>
  )
}

export default PairpotFooter
