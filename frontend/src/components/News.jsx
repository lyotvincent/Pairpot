import 'animate.css'
import contactImg from "../assets/img/contact.png"
import { ConfigProvider, Divider, Layout, Space } from 'antd'
import { Card, Row, Col } from 'antd'
import { List } from 'echarts'
const { Content, Footer } = Layout
import githubLogo from "../assets/img/github.svg"

const News = () => {
    return (
        <ConfigProvider theme={{
          components: {
            Card: {
              colorBgContainer: "#fff",
              headerFontSize: 32,
              headerHeight: 48,
              borderRadiusLG: 12,
              colorTextDescription: "rgba(0, 0, 0, 0.6)",
            }
          }
        }}>
          <Layout>
            <Content>
              <section className="banner" id="help">
                <Row>
                    <Col>
                        <h1>Latest News</h1>
                    </Col>
                    <Divider/>

                  <Col >
                    <div>
                      <h2>
                        Citation
                      </h2>
                      <h4 style={{fontFamily:'initial'}}>
                        Ruan Z., Lin F., Zhang Z., et al. (2024) Pairpot: a database with real-time lasso-based analysis tailored for paired single-cell and spatial transcriptomics. Nucleic Acids Research, gkae986, &nbsp;
                        <a href="https://doi.org/10.1093/nar/gkae986"target='_blank'
                            style={{
                            color: 'rgba(0, 0, 0, 0.6)',
                            textDecoration: 'none',
                            fontSize: 24
                            }}>
                            https://doi.org/10.1093/nar/gkae986
                        </a>
                      </h4>
                    </div>
                  </Col>
                
                  <Divider/>

                  <Col>
                    <div>
                        <h2>
                            Web Updating News
                        </h2>
                        <h4>
                            {/* <span style={{color:'gray'}}>2024.11.7&nbsp;&nbsp;</span>
                            We have fixed the calculation issue of cell type proportion in Pairview for some datasets.
                            <br/> */}
                            <span style={{color:'gray'}}>2024.11.4&nbsp;&nbsp;</span>
                            Pairpot online formally. 
                        </h4>
                    </div>
                  </Col>
                </Row>
                
              </section>
            </Content>
    
            <Footer>
              <Row>
    
              </Row>
              <Row>
                <Col span={4}>
                  <Space direction='horizontal'>
                    <h3>Contact us
                      {' '}
                      <img src={contactImg} style={{ width: '20%', height: '20%' }} alt="contact"></img>
                    </h3>
                  </Space>
                  <a href="https://beian.miit.gov.cn" target='_blank'
                    style={{
                      color: 'rgba(0, 0, 0, 0.6)',
                      textDecoration: 'none',
                      fontSize: 13
                    }}>
                    津ICP备2024022781号-1
                  </a>
                </Col>
                <Col>
                <a href="https://github.com/lyotvincent/Pairpot">
                    <img src={githubLogo} alt="GitHub Logo" />
                </a>
                </Col>
                <Col offset={1}>
                  <p>
                    *Zhihan Ruan, Centre for Bioinformatics and Intelligent
                    Medicine, Nankai University, rrrzhan@nankai.edu.cn
                  </p>
                  <p>
                    *Jian Liu, State Key Laboratory of Medical Chemical Biology, College of Computer Science, Nankai University, jianliu@nankai.edu.cn
                  </p>
                </Col>
              </Row>
            </Footer>
          </Layout>
        </ConfigProvider>
      )
}
export default News