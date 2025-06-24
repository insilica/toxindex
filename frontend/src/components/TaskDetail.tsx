import React, { useEffect, useState, useCallback } from "react";
import { useParams, useNavigate } from "react-router-dom";
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';

interface Task {
  task_id: string;
  title: string;
  archived?: boolean;
  environment_id?: string;
  created_at?: string;
  finished_at?: string;
  session_id?: string;
  user_id?: string;
  workflow_id?: number;
}

const getDuration = (start?: string, end?: string) => {
  if (!start || !end) return null;
  const startTime = new Date(start).getTime();
  const endTime = new Date(end).getTime();
  const diff = Math.max(0, endTime - startTime);
  const mins = Math.floor(diff / 60000);
  const secs = Math.floor((diff % 60000) / 1000);
  return `${mins}m ${secs}s`;
};

const TaskDetail: React.FC = () => {
  const { task_id } = useParams<{ task_id: string }>();
  const [task, setTask] = useState<Task | null>(null);
  const [loading, setLoading] = useState(true);
  const [assistantMessage, setAssistantMessage] = useState<string | null>(null);
  const [messageLoading, setMessageLoading] = useState(false);
  const navigate = useNavigate();

  const fetchTask = useCallback(async () => {
    if (!task_id) return;
    setLoading(true);
    try {
      const res = await fetch(`/api/tasks/${task_id}`, { credentials: "include" });
      if (!res.ok) {
        setTask(null);
        setLoading(false);
        return;
      }
      const data = await res.json();
      if (data.error) {
        setTask(null);
      } else {
        setTask(data);
      }
    } catch (e) {
      setTask(null);
    } finally {
      setLoading(false);
    }
  }, [task_id]);

  // Fetch assistant message for this task
  useEffect(() => {
    const fetchAssistantMessage = async () => {
      if (!task_id) return;
      setMessageLoading(true);
      try {
        const res = await fetch(`/api/tasks/${task_id}/messages`, { credentials: "include" });
        if (!res.ok) throw new Error('Failed to fetch messages');
        const data = await res.json();
        if (data.messages && Array.isArray(data.messages)) {
          const assistantMsg = data.messages.find((m: any) => m.role === "assistant");
          if (assistantMsg) {
            setAssistantMessage(assistantMsg.content);
            return;
          }
        }
        setAssistantMessage(null); // No assistant message found
      } catch (e) {
        setAssistantMessage('Failed to load assistant message.');
        if (e instanceof Error) {
          console.error('Error loading assistant message:', e.message);
        } else {
          console.error('Unknown error loading assistant message:', e);
        }
      } finally {
        setMessageLoading(false);
      }
    };
    fetchAssistantMessage();
  }, [task_id]);

  useEffect(() => {
    fetchTask();
  }, [fetchTask]);

  if (loading) return <div className="text-white p-8">Loading...</div>;
  if (!task) return <div className="text-white p-8">Task not found.</div>;

  const duration = getDuration(task.created_at, task.finished_at);

  return (
    <div className="max-w-7xl mx-auto p-8 bg-gray-900 rounded-lg shadow text-white mt-12 flex flex-col min-h-[70vh]">
      <h2 className="text-3xl font-bold mb-6 text-center">{task.title}</h2>
      <div className="flex-1 flex flex-col items-center justify-center">
        {messageLoading ? (
          <div className="text-gray-400">Loading result...</div>
        ) : assistantMessage ? (
          <div className="prose prose-invert w-full max-w-full bg-gray-800 rounded-lg p-6 shadow-inner overflow-x-auto" style={{ lineHeight: 2 }}>
            <ReactMarkdown remarkPlugins={[remarkGfm]}>{assistantMessage}</ReactMarkdown>
          </div>
        ) : (
          <div className="text-gray-400">No assistant message found for this task.</div>
        )}
      </div>
      <footer className="mt-8 pt-6 border-t border-gray-700 grid grid-cols-1 md:grid-cols-2 gap-4 text-sm bg-gray-950 rounded-b-lg">
        <div>
          <b>Tool:</b> {task.workflow_id === 1 ? "ToxIndex RAP" : task.workflow_id === 2 ? "ToxIndex Vanilla" : task.workflow_id ?? <span className="text-gray-400">Unknown</span>}
        </div>
        <div>
          <b>Executed by user:</b> {task.user_id ? (
            <button
              className="text-purple-400 underline hover:text-purple-300 focus:outline-none bg-transparent p-0 shadow-none border-none"
              style={{ background: 'none', boxShadow: 'none', border: 'none' }}
              onClick={() => navigate(`/user/${task.user_id}`)}
              title="Go to user profile"
              aria-label="Go to user profile"
            >
              {task.user_id}
            </button>
          ) : (
            <span className="text-gray-400">Unknown</span>
          )}
        </div>
        <div>
          <b>Created:</b> {task.created_at ? new Date(task.created_at).toLocaleString() : <span className="text-gray-400">Unknown</span>}
        </div>
        <div>
          <b>Finished:</b> {task.finished_at ? new Date(task.finished_at).toLocaleString() : <span className="text-gray-400">-</span>}
        </div>
        <div>
          <b>Duration:</b> {duration ?? <span className="text-gray-400">-</span>}
        </div>
        <div>
          <b>Environment:</b> {task.environment_id ? (
            <button
              className="text-blue-400 underline hover:text-blue-300 focus:outline-none bg-transparent p-0 shadow-none border-none"
              style={{ background: 'none', boxShadow: 'none', border: 'none' }}
              onClick={() => navigate(`/environment/${task.environment_id}`)}
              title="Go to environment"
              aria-label="Go to environment"
            >
              {task.environment_id}
            </button>
          ) : (
            <span className="text-gray-400">None</span>
          )}
        </div>
        <div>
          <b>Chat Session:</b> {task.session_id ? (
            <button
              className="text-green-400 underline hover:text-green-300 focus:outline-none bg-transparent p-0 shadow-none border-none"
              style={{ background: 'none', boxShadow: 'none', border: 'none' }}
              onClick={() => navigate(`/chat/session/${task.session_id}`)}
              title="Go to chat session"
              aria-label="Go to chat session"
            >
              {task.session_id}
            </button>
          ) : (
            <span className="text-gray-400">None</span>
          )}
        </div>
        <div className="col-span-1 md:col-span-2 flex justify-end">
          <button className="mt-2 px-4 py-2 bg-green-700 rounded" onClick={() => navigate(-1)} aria-label="Back">Back</button>
        </div>
      </footer>
    </div>
  );
};

export default TaskDetail; 