import React, { useState } from 'react'
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
  DatePicker,
  TimePicker,
} from 'antd'

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

const { Option } = Select
const formItemLayout = {
  labelCol: {
    span: 4,
  },
  wrapperCol: {
    span: 18,
  },
}

const SubmitMeta = (props) => {
  const { onFinish } = props
  const [timeType1, setTimeType1] = useState('time')
  const [timeType2, setTimeType2] = useState('time')
  const PickerWithType = ({ type, onChange }) => {
    if (type === 'time') return <TimePicker onChange={onChange} />
    if (type === 'date') return <DatePicker onChange={onChange} />
    return <DatePicker picker={type} onChange={onChange} />
  }

  return (
    <>
      <h2 style={{ color: 'black' }}>Submit New Dataset Metadata</h2>
      <br />
      <Form
        name="ST reference dataset"
        {...formItemLayout}
        onFinish={onFinish}
        initialValues={initialValues}
        style={{
          maxWidth: 1200,
        }}>
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
            maxLength={500}
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
          <Select placeholder="Please select the technologies" mode="multiple">
            {technologiesOption.map((item) => (
              <Option value={item}>{item}</Option>
            ))}
          </Select>
        </Form.Item>
        <Form.Item {...formItemLayout} name="organ_parts" label="Organ Parts">
          <Input placeholder="Please input the Organ Parts if available." />
        </Form.Item>
        {/*disease*/}
        <Form.Item {...formItemLayout} name="disease" label="Disease">
          <Input placeholder="Please input the disease" />
        </Form.Item>
        {/*development_stages*/}
        <Form.Item
          {...formItemLayout}
          name="development_stages"
          label="Development stages">
          <Input placeholder="Please input the development stages if available." />
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
              labelCol={{ span: 8 }}
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
              labelCol={{ span: 8 }}
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
        <Row>
          <Col span={12}>
            <Form.Item
              labelCol={{ span: 8 }}
              wrapperCol={{
                span: 12,
                offset: 0,
              }}
              name="submission_date"
              label="Submission Date">
              <Space>
                <Select value={timeType1} onChange={setTimeType1}>
                  <Option value="time">Time</Option>
                  <Option value="date">Date</Option>
                  <Option value="week">Week</Option>
                  <Option value="month">Month</Option>
                  <Option value="quarter">Quarter</Option>
                  <Option value="year">Year</Option>
                </Select>
                <PickerWithType
                  type={timeType1}
                  onChange={(value) => console.log(value)}
                />
              </Space>
            </Form.Item>
          </Col>
          <Col span={12}>
            <Form.Item
              labelCol={{ span: 6 }}
              wrapperCol={{
                span: 12,
                offset: 0,
              }}
              name="last_modified"
              label="Last Modified">
              <Space>
                <Select value={timeType2} onChange={setTimeType2}>
                  <Option value="time">Time</Option>
                  <Option value="date">Date</Option>
                  <Option value="week">Week</Option>
                  <Option value="month">Month</Option>
                  <Option value="quarter">Quarter</Option>
                  <Option value="year">Year</Option>
                </Select>
                <PickerWithType
                  type={timeType2}
                  onChange={(value) => console.log(value)}
                />
              </Space>
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
        <Form.Item {...formItemLayout} name="pmid" label="PMID">
          <Input placeholder="Please input the PMID" />
        </Form.Item>
        <Form.Item {...formItemLayout} name="platforms" label="Platforms">
          <Input placeholder="Please input the platforms" />
        </Form.Item>

        <Form.Item
          wrapperCol={{
            span: 12,
            offset: 6,
          }}>
          <Space>
            <Button htmlType="submit" type="primary">
              Submit
            </Button>
            <Button htmlType="reset">reset</Button>
          </Space>
        </Form.Item>
      </Form>
    </>
  )
}

export default SubmitMeta
