import React, { useEffect, useState } from 'react'
import {
  Button,
  Input,
  Col,
  Form,
  InputNumber,
  Radio,
  Row,
  Select,
  Space,
} from 'antd'
import { useNavigate } from 'react-router-dom'

const initialValues = {
  dataset_id: '',
  title: '',
  contributors: '',
  summary: '',
  overall_design: '',
  species: [],
  tissues: [],
  technologies: [],
  disease: '',
  sex: 'Unknown',
  n_samples: 0,
  cells: 0,
  spots: 0,
  genes: 0,
  contacts: '',
  citation: '',
  accessions: '',
  platforms: '',
  PMID: '',
}

const speciesOption = [
  'Ambystoma mexicanum',
  'Arabidopsis thaliana',
  'Caenorhabditis elegans',
  'Canis lupus familiaris',
  'Danio rerio',
  'Drosophila melanogaster',
  'Gallus gallus',
  'Glycine max',
  'Homo sapiens',
  'Macaca fascicularis',
  'Maize',
  'Mus musculus',
  'Oryctolagus cuniculus',
  'Phalaenopsis aphrodite',
  'Rattus norvegicus',
  'Sus scrofa',
  'Xenopus laevis',
  'Xenopus tropicalis',
  'Unknown',
]

const tissuesOption = [
  'Basal cell',
  'Bladder',
  'Blood',
  'Blood vessel',
  'Bone',
  'Bone marrow',
  'Brain',
  'Breast',
  'Cell line',
  'Cerebella',
  'Cervical',
  'Colon',
  'Cortex',
  'Dorsal aorta',
  'Ear',
  'Embryo',
  'Epithelium',
  'Eye',
  'Eymbro',
  'Flower',
  'Gastrula',
  'Gland',
  'Heart',
  'Ileum',
  'Intestine',
  'Kidney',
  'Lacrimal Gland',
  'Leaf',
  'Legs',
  'Liver',
  'Lung',
  'Lymph',
  'Mucosa',
  'Muscle',
  'Neuron',
  'Ovarian',
  'Pancreas',
  'Pituitary',
  'Plant part',
  'Prostate',
  'Skeleton',
  'Skin',
  'Spinal cord',
  'Spleen',
  'Stomach',
  'Synovial',
  'Thymus',
  'Tumor',
  'Uterus',
  'Unknown',
]

const technologiesOption = [
  '10x Visium',
  'APEX-seq',
  'CUT&Tag',
  'ClampFISH',
  'ClumpSeq',
  'DBiT-seq',
  'EASI-FISH',
  'Geo-seq',
  'GeoMx DSP',
  'HDST',
  'LCM-seq',
  'MERFISH',
  'PIC-seq',
  'Proximity RNA-seq',
  'RNAseq',
  'STARmap',
  'Seq-Scope',
  'Slide-seqV2',
  'Spatial Transcriptomics',
  'Stereo-Seq',
  'Tomo-seq',
  'XYZeq',
  'scATAC',
  'scRNA',
  'sci-ATAC-seq',
  'sci-Space',
  'sciMAP-ATAC-seq',
  'seqFISH+',
  'seqFish',
  'snRNA',
]

const sexOption = ['Male', 'Female', 'Plant', 'Other', 'Unknown']

const pairOption = [
  'From original dataset',
  'From existing dataset',
  'From external dataset',
  'From Atlas',
]

const rateOption = []

const { Option } = Select
const formItemLayout = {
  labelCol: {
    span: 5,
  },
  wrapperCol: {
    span: 18,
  },
}

const SubmitLink = (props) => {
  const { onFinish, location, navigate, submitForm } = props
  const HandlePairState = (value) => {
    let fieldsValue = {}
    if (value === 'From original dataset') {
      fieldsValue = location.state
      fieldsValue['has_paired'] = fieldsValue.dataset_id
    } else if (value === 'From external dataset') {
      fieldsValue = initialValues
      fieldsValue['has_paired'] = location.state.dataset_id
    }
    submitForm.setFieldsValue(fieldsValue)
  }
  return (
    <>
      <h2 style={{ color: 'black' }}>Link Paired Dataset Metadata</h2>
      <br />
      <Row gutter={[0, 0]}>
        {/*left is reference form existing ST metadata*/}
        <Col flex={2}>
          <h4 style={{ color: 'black', textAlign: 'center' }}>
            Reference Dataset
          </h4>
          <Form
            name="ST reference dataset"
            {...formItemLayout}
            onFinish={onFinish}
            initialValues={location.state}
            style={{
              maxWidth: 1200,
            }}>
            <Form.Item label="Dataset ID" name="dataset_id" required={true}>
              <Input placeholder="Please input the dataset_id" />
            </Form.Item>
            <Form.Item
              {...formItemLayout}
              name="title"
              label="Title"
              rules={[
                {
                  required: true,
                  message: 'The title is required.',
                },
              ]}>
              <Input.TextArea
                placeholder="Please input the title"
                autoSize={{ minRows: 1, maxRows: 3 }}
              />
            </Form.Item>
            <Form.Item
              {...formItemLayout}
              name="contributors"
              label="Contributors"
              required>
              <Input.TextArea
                placeholder="Please add the contributors"
                autoSize={{ minRows: 1, maxRows: 3 }}
              />
            </Form.Item>
            <Form.Item name="summary" label="Summary" required>
              <Input.TextArea
                showCount
                maxLength={5000}
                placeholder="Please enter the summary"
                autoSize={{ minRows: 3, maxRows: 6 }}
              />
            </Form.Item>
            <Form.Item name="overall_design" label="Overall Design" required>
              <Input.TextArea
                showCount
                maxLength={5000}
                placeholder="Please enter the overall design"
                autoSize={{ minRows: 3, maxRows: 6 }}
              />
            </Form.Item>
            {/* species */}
            <Form.Item
              name="species"
              label="Species"
              hasFeedback
              rules={[
                {
                  required: true,
                  message: 'The species are required.',
                },
              ]}>
              <Select placeholder="Please select the species" mode="multiple">
                {speciesOption.map((item) => (
                  <Option value={item}>{item}</Option>
                ))}
              </Select>
            </Form.Item>
            {/* tissues */}
            <Form.Item
              name="tissues"
              label="Tissues"
              hasFeedback
              rules={[
                {
                  required: true,
                  message: 'The Tissues are required.',
                },
              ]}>
              <Select placeholder="Please select the tissues" mode="multiple">
                {tissuesOption.map((item) => (
                  <Option value={item}>{item}</Option>
                ))}
              </Select>
            </Form.Item>
            {/* technologies */}
            <Form.Item
              name="technologies"
              label="Technologies"
              hasFeedback
              rules={[
                {
                  required: true,
                  message: 'The Technologies are required.',
                },
              ]}>
              <Select
                placeholder="Please select the technologies"
                mode="multiple">
                {technologiesOption.map((item) => (
                  <Option value={item}>{item}</Option>
                ))}
              </Select>
            </Form.Item>
            {/* sex */}
            <Form.Item
              name="sex"
              label="Sex"
              hasFeedback
              rules={[
                {
                  required: true,
                  message: 'The Sex is required.',
                },
              ]}>
              <Radio.Group>
                {sexOption.map((item) => (
                  <Radio value={item}>{item}</Radio>
                ))}
              </Radio.Group>
            </Form.Item>
            <Row>
              <Col span={12}>
                <Form.Item
                  labelCol={{ span: 12 }}
                  wrapperCol={{
                    span: 6,
                    offset: 0,
                  }}
                  label="Samples"
                  name="n_samples"
                  rules={[
                    {
                      required: true,
                      message: 'The number of samples are required.',
                    },
                  ]}>
                  <InputNumber />
                </Form.Item>
              </Col>
              <Col span={12}>
                <Form.Item
                  labelCol={{ span: 6 }}
                  wrapperCol={{
                    span: 6,
                    offset: 0,
                  }}
                  label="Cells"
                  name="cells">
                  <InputNumber />
                </Form.Item>
              </Col>
            </Row>
            <Row>
              <Col span={12}>
                <Form.Item
                  labelCol={{ span: 12 }}
                  wrapperCol={{
                    span: 6,
                  }}
                  label="Spots"
                  name="spots">
                  <InputNumber />
                </Form.Item>
              </Col>
              <Col span={12}>
                <Form.Item
                  labelCol={{ span: 6 }}
                  wrapperCol={{
                    span: 6,
                    offset: 0,
                  }}
                  label="Genes"
                  name="genes">
                  <InputNumber />
                </Form.Item>
              </Col>
            </Row>

            <Form.Item {...formItemLayout} name="contacts" label="Contacts">
              <Input placeholder="Please input the contacts" />
            </Form.Item>

            <Form.Item {...formItemLayout} name="citation" label="Citation">
              <Input placeholder="Please input the citation" />
            </Form.Item>
            <Form.Item
              {...formItemLayout}
              name="accessions"
              label="Accessions"
              hasFeedback
              rules={[
                {
                  required: true,
                  message: 'The accession URLs are required.',
                },
              ]}>
              <Input type="url" placeholder="Please input the accession URLs" />
            </Form.Item>
            <Form.Item {...formItemLayout} name="platforms" label="Platforms">
              <Input placeholder="Please input the platforms" />
            </Form.Item>
            <Form.Item {...formItemLayout} name="PMID" label="PMID">
              <Input placeholder="Please input the PMID" />
            </Form.Item>

            <Form.Item
              wrapperCol={{
                span: 12,
                offset: 6,
              }}>
              <Space>
                <Button
                  type="primary"
                  onClick={() => {
                    navigate('/database')
                  }}>
                  Find a ST dataset
                </Button>
                <Button>Correct Meta</Button>
                <Button htmlType="reset">reset</Button>
              </Space>
            </Form.Item>
          </Form>
        </Col>
        <Col flex={2}>
          <h4 style={{ color: 'black', textAlign: 'center' }}>Link Dataset</h4>
          <Form
            name="submit-paired-dataset"
            {...formItemLayout}
            onFinish={onFinish}
            form={submitForm}
            style={{
              maxWidth: 1200,
            }}>
            <Form.Item
              {...formItemLayout}
              name="pair-state"
              label="Pair State"
              rules={[
                {
                  required: true,
                  message: 'Pair state is required.',
                },
              ]}>
              <Select
                placeholder="Please select the pair states"
                onChange={HandlePairState}>
                {pairOption.map((item) => (
                  <Option value={item}>{item}</Option>
                ))}
              </Select>
            </Form.Item>
            <Form.Item label="Paired ID" name="has_paired" required={true}>
              <Input placeholder="Please input the paired dataset_id" />
            </Form.Item>
            <Form.Item
              {...formItemLayout}
              name="title"
              label="Title"
              rules={[
                {
                  required: true,
                  message: 'The title is required.',
                },
              ]}>
              <Input.TextArea
                placeholder="Please input the title"
                autoSize={{ minRows: 1, maxRows: 3 }}
              />
            </Form.Item>
            <Form.Item
              {...formItemLayout}
              name="contributors"
              label="Contributors"
              required>
              <Input placeholder="Please add the contributors" />
            </Form.Item>
            <Form.Item name="summary" label="Summary" required>
              <Input.TextArea
                showCount
                maxLength={5000}
                placeholder="Please enter the summary"
              />
            </Form.Item>
            <Form.Item name="overall_design" label="Overall Design" required>
              <Input.TextArea
                showCount
                maxLength={5000}
                placeholder="Please enter the overall design"
              />
            </Form.Item>
            {/* species */}
            <Form.Item
              name="species"
              label="Species"
              hasFeedback
              rules={[
                {
                  required: true,
                  message: 'The species are required.',
                },
              ]}>
              <Select placeholder="Please select the species" mode="multiple">
                {speciesOption.map((item) => (
                  <Option value={item}>{item}</Option>
                ))}
              </Select>
            </Form.Item>
            {/* tissues */}
            <Form.Item
              name="tissues"
              label="Tissues"
              hasFeedback
              rules={[
                {
                  required: true,
                  message: 'The Tissues are required.',
                },
              ]}>
              <Select placeholder="Please select the tissues" mode="multiple">
                {tissuesOption.map((item) => (
                  <Option value={item}>{item}</Option>
                ))}
              </Select>
            </Form.Item>
            {/* technologies */}
            <Form.Item
              name="technologies"
              label="Technologies"
              hasFeedback
              rules={[
                {
                  required: true,
                  message: 'The Technologies are required.',
                },
              ]}>
              <Select
                placeholder="Please select the technologies"
                mode="multiple">
                {technologiesOption.map((item) => (
                  <Option value={item}>{item}</Option>
                ))}
              </Select>
            </Form.Item>
            {/* sex */}
            <Form.Item
              name="sex"
              label="Sex"
              hasFeedback
              rules={[
                {
                  required: true,
                  message: 'The Sex is required.',
                },
              ]}>
              <Radio.Group>
                {sexOption.map((item) => (
                  <Radio value={item}>{item}</Radio>
                ))}
              </Radio.Group>
            </Form.Item>
            <Row>
              <Col span={12}>
                <Form.Item
                  labelCol={{ span: 10 }}
                  wrapperCol={{
                    span: 8,
                    offset: 0,
                  }}
                  label="Samples"
                  name="n_samples"
                  rules={[
                    {
                      required: true,
                      message: 'The number of samples are required.',
                    },
                  ]}>
                  <InputNumber />
                </Form.Item>
              </Col>
              <Col span={12}>
                <Form.Item
                  labelCol={{ span: 5 }}
                  wrapperCol={{
                    span: 8,
                    offset: 0,
                  }}
                  label="Cells"
                  name="cells">
                  <InputNumber />
                </Form.Item>
              </Col>
            </Row>
            <Row>
              <Col span={12}>
                <Form.Item
                  labelCol={{ span: 10 }}
                  wrapperCol={{
                    span: 8,
                    offset: 0,
                  }}
                  label="Spots"
                  name="spots">
                  <InputNumber />
                </Form.Item>
              </Col>
              <Col span={12}>
                <Form.Item
                  labelCol={{ span: 5 }}
                  wrapperCol={{
                    span: 8,
                    offset: 0,
                  }}
                  label="Genes"
                  name="genes">
                  <InputNumber />
                </Form.Item>
              </Col>
            </Row>

            <Form.Item {...formItemLayout} name="contacts" label="Contacts">
              <Input placeholder="Please input the contacts" />
            </Form.Item>

            <Form.Item {...formItemLayout} name="citation" label="Citation">
              <Input placeholder="Please input the citation" />
            </Form.Item>
            <Form.Item
              {...formItemLayout}
              name="accessions"
              label="Accessions"
              hasFeedback
              rules={[
                {
                  required: true,
                  message: 'The accession URLs are required.',
                },
              ]}>
              <Input.TextArea
                autoSize={{ minRows: 1, maxRows: 3 }}
                placeholder="Please input the accession URLs"
              />
            </Form.Item>
            <Form.Item {...formItemLayout} name="platforms" label="Platforms">
              <Input placeholder="Please input the platforms" />
            </Form.Item>

            <Form.Item {...formItemLayout} name="PMID" label="PMID">
              <Input placeholder="Please input the PMID" />
            </Form.Item>
            
            <Form.Item
              wrapperCol={{
                span: 12,
                offset: 6,
              }}>
              <Space>
                <Button type="primary" htmlType="submit">
                  Submit
                </Button>
                <Button htmlType="reset">reset</Button>
              </Space>
            </Form.Item>
          </Form>
        </Col>
      </Row>
    </>
  )
}

export default SubmitLink
