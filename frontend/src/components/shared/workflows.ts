// List of workflows

export const WORKFLOWS = [
  { id: 'toxindex-rap', label: 'ToxIndex RAP' },
  { id: 'toxindex-vanilla', label: 'ToxIndex Vanilla' },
  { id: 'toxindex-pathway', label: 'ToxIndex Pathway' },
  { id: 'toxindex-ChemParser', label: 'ChemParser' },
  { id: 'toxindex-5th', label: 'ToxIndex 5th' }, // unsupported
];

// Only active workflows for id mapping
const ACTIVE_WORKFLOWS = WORKFLOWS.filter(w => w.id !== 'toxindex-5th');

export function getWorkflowId(selectedModel: string): number {
  const idx = ACTIVE_WORKFLOWS.findIndex(w => w.id === selectedModel);
  return idx === -1 ? 0 : idx + 1; // 1-based index, 0 if not found/unsupported
}

export function getWorkflowLabelById(workflow_id: number): string {
  // workflow_id is 1-based index into ACTIVE_WORKFLOWS
  if (workflow_id > 0 && workflow_id <= ACTIVE_WORKFLOWS.length) {
    return ACTIVE_WORKFLOWS[workflow_id - 1].label;
  }
  return 'Unknown';
} 