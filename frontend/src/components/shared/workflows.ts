// Dynamically loads workflows from resources/default_workflows.json
// This keeps frontend and backend in sync automatically

export interface Workflow {
  workflow_id: number;
  frontend_id: string;
  title: string;
  label: string;
  description: string;
  initial_prompt: string;
  celery_task: string;
}

// Cache for loaded workflows
let workflowsCache: Workflow[] | null = null;
let isLoading = false;
let loadPromise: Promise<Workflow[]> | null = null;

// Load workflows from JSON file
export async function loadWorkflows(forceRefresh = false): Promise<Workflow[]> {
  if (workflowsCache && !forceRefresh) {
    return workflowsCache;
  }

  if (loadPromise) {
    return loadPromise;
  }

  if (isLoading) {
    // Wait for the current load to complete
    while (isLoading) {
      await new Promise(resolve => setTimeout(resolve, 10));
    }
    return workflowsCache || [];
  }

  isLoading = true;
  loadPromise = (async () => {
    try {
      // Read from resources directory via API endpoint
      const response = await fetch('/api/workflows/config');
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      const data = await response.json();
      workflowsCache = data.workflows as Workflow[];
      return workflowsCache;
    } catch (error) {
      console.error('Failed to load workflows from JSON:', error);
      // Return empty array as fallback
      return [];
    } finally {
      isLoading = false;
      loadPromise = null;
    }
  })();

  return loadPromise;
}

// Get workflows (synchronous wrapper that uses cached data)
export function getWorkflows(): Workflow[] {
  return workflowsCache || [];
}

// Clear workflows cache to force reload
export function clearWorkflowsCache(): void {
  workflowsCache = null;
  isLoading = false;
  loadPromise = null;
}

// Initialize workflows on module load
export async function initializeWorkflows(): Promise<void> {
  await loadWorkflows();
}

// Legacy compatibility - returns current workflows (will be empty initially)
export const WORKFLOWS = getWorkflows();

export function getWorkflowId(selectedModel: string): number {
  const workflows = getWorkflows();
  const workflow = workflows.find(w => w.frontend_id === selectedModel);
  return workflow ? workflow.workflow_id : 0;
}

export function getWorkflowLabelById(workflow_id: number): string {
  const workflows = getWorkflows();
  const workflow = workflows.find(w => w.workflow_id === workflow_id);
  return workflow ? workflow.label : 'Unknown';
}

export function getWorkflowById(workflow_id: number): Workflow | undefined {
  const workflows = getWorkflows();
  return workflows.find(w => w.workflow_id === workflow_id);
}

export function getWorkflowByFrontendId(frontend_id: string): Workflow | undefined {
  const workflows = getWorkflows();
  return workflows.find(w => w.frontend_id === frontend_id);
}

// Initialize workflows when this module is imported
// Removed automatic initialization to prevent stale cache 